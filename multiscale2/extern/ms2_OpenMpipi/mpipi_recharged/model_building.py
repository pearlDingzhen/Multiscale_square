"""
Model building utilities for OpenMpipi.

Provides functions for:
- Calculating target box vectors
- Building combined models from multiple chains
- Equilibrating slab configurations
"""

import numpy as np
import openmm.app as app
import openmm as mm
import openmm.unit as unit
import scipy.constants as c
import os
import random
import warnings

from .system_building import get_mpipi_system
from .constants import PLATFORM, PROPERTIES


# ====== CALVADOS-style placement functions ======

def draw_starting_vec(box):
    """
    Draw random position within the simulation box.
    
    Args:
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.
    
    Returns:
        numpy.ndarray: Random 3D position within the box.
    """
    return np.random.random(3) * np.array(box)


def check_walls(x, box):
    """
    Check if molecule position is within the simulation box.
    
    Args:
        x (numpy.ndarray): Positions of all atoms in the molecule.
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.
    
    Returns:
        bool: True if molecule is outside the box (clashes with wall), False otherwise.
    """
    box = np.array(box)
    if np.min(x) < 0:
        return True  # clash with left box wall
    d = box - x
    if np.min(d) < 0:
        return True  # clash with right box wall
    return False  # molecule in box


def check_clash(x, pos_others, box, cutoff=0.75):
    """
    Check for clashes with other particles using distance array.
    
    Args:
        x (numpy.ndarray): Positions of atoms in the new molecule.
        pos_others (list): List of position arrays for other molecules already placed.
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.
        cutoff (float): Minimum distance (in nm) to avoid clash. Default 0.75 nm.
    
    Returns:
        bool: True if clash detected, False otherwise.
    """
    boxfull = np.append(box, [90, 90, 90])  # Add angles for distance calculation
    
    if len(pos_others) == 0:
        return False  # no other particles
    
    # Convert to numpy array
    pos_others_arr = np.array(pos_others)
    if len(pos_others_arr) == 0:
        return False
    
    # Calculate minimum distances
    min_dists = _distance_array(x, pos_others_arr, boxfull)
    
    # min_dists is (K, N) where K is number of other molecules, N is atoms in new molecule
    # We check if ANY atom in the new molecule is too close to ANY other molecule
    if min_dists.ndim == 2:
        # Check if any molecule has all its atoms far enough away
        # Actually we want: is ANY atom in x too close to ANY atom in any other molecule?
        # So we check if min distance for any atom to any molecule is < cutoff
        overall_min = np.min(min_dists)
    else:
        overall_min = np.min(min_dists)
    
    if overall_min < cutoff:
        return True  # clash with other particles
    else:
        return False  # no clash


def _distance_array(x, y, box):
    """
    Calculate minimum distance from each point in x to any point in y.
    
    Args:
        x (numpy.ndarray): Coordinates of first set (N, 3).
        y (numpy.ndarray): Coordinates of second set (M, 3) or (K, M, 3) for multiple molecules.
        box (numpy.ndarray): Box dimensions [Lx, Ly, Lz, alpha, beta, gamma].
    
    Returns:
        numpy.ndarray: Minimum distance array (N,) if y is 2D, or (K, N) if y is 3D.
    """
    box = np.array(box)
    Lx, Ly, Lz = box[:3]
    box_arr = np.array([Lx, Ly, Lz])
    
    # Handle y shape: if 3D (multiple molecules), flatten to 2D for distance calc
    if y.ndim == 3:
        # y shape: (K, M, 3) - K molecules, M atoms each
        K, M, _ = y.shape
        y_flat = y.reshape(-1, 3)  # (K*M, 3)
        
        # Calculate distances for each molecule separately
        result = np.zeros((K, x.shape[0]))
        for k in range(K):
            y_k = y[k]  # (M, 3)
            dx = x[:, np.newaxis, :] - y_k[np.newaxis, :, :]  # (N, M, 3)
            dx = dx - box_arr[np.newaxis, np.newaxis, :] * np.round(dx / box_arr[np.newaxis, np.newaxis, :])
            dist = np.sqrt(np.sum(dx**2, axis=2))  # (N, M)
            result[k] = np.min(dist, axis=1)  # (N,) - min distance for each atom in x
        return result  # (K, N)
    else:
        # y shape: (M, 3) - single molecule
        dx = x[:, np.newaxis, :] - y[np.newaxis, :, :]  # (N, M, 3)
        dx = dx - box_arr[np.newaxis, np.newaxis, :] * np.round(dx / box_arr[np.newaxis, np.newaxis, :])
        dist = np.sqrt(np.sum(dx**2, axis=2))  # (N, M)
        return np.min(dist, axis=1)  # (N,)


def random_placement(box, pos_others, xinit, ntries=10000, cutoff=0.75, verbose=False):
    """
    Randomly place a molecule within the box, avoiding clashes and walls.
    
    Args:
        box (array-like): Box dimensions [Lx, Ly, Lz] in nm.
        pos_others (list): List of position arrays for other molecules already placed.
        xinit (numpy.ndarray): Initial coordinates of the molecule to place (centered at origin).
        ntries (int): Maximum number of placement attempts. Default 10000.
        cutoff (float): Minimum distance to avoid clash. Default 0.75 nm.
        verbose (bool): Print verbose output.
    
    Returns:
        numpy.ndarray: Final positions of the placed molecule.
    
    Raises:
        ValueError: If unable to place molecule within ntries attempts.
    """
    for ntry in range(ntries):
        x0 = draw_starting_vec(box)  # random point in box
        xs = x0 + xinit  # shift molecule to that position
        
        # Check if within box walls
        if check_walls(xs, box):
            if verbose and ntry % 1000 == 0:
                print(f'  Try {ntry}: outside box walls', flush=True)
            continue
        
        # Check for clashes with other molecules
        if check_clash(xs, pos_others, box, cutoff=cutoff):
            if verbose and ntry % 1000 == 0:
                print(f'  Try {ntry}: clash detected', flush=True)
            continue
        
        # Success!
        if verbose and ntry > 0:
            print(f'  Placed after {ntry + 1} attempts', flush=True)
        return xs
    
    raise ValueError(f'Failed to place molecule after {ntries} attempts. '
                     f'Consider increasing box size or reducing molecule count.')


# ====== Box vector calculation ======

def calculate_target_box_vectors(chain_info, 
                                 box_size=None,
                                 target_density=0.1*unit.gram/unit.centimeter**3):
    """
    Calculate the initial target box vectors for the system.
    
    If box_size is provided, use it directly.
    Otherwise, calculate from total mass and target density.

    Args:
        chain_info (dict): 
            A dictionary with keys = chain objects and values = number of copies to place.
            Each chain object must have:
              - chain_mass (in daltons)
              - min_rg_coords (initial coordinates)
        box_size (list, optional): 
            Box size [Lx, Ly, Lz] in nm. If provided, use this directly.
        target_density (Quantity, optional): 
            Target density as an OpenMM `Quantity` (default 0.1 g/cm^3).

    Returns:
        Quantity: 
            A 3x3 array (in nm) of box vectors suitable for periodic boundary conditions.
    """
    if box_size is not None:
        # Use provided box size directly
        Lx, Ly, Lz = box_size
        return np.array([[Lx, 0, 0], 
                         [0, Ly, 0], 
                         [0, 0, Lz]]) * unit.nanometer
    
    # Calculate from total mass and target density
    total_mass_g = sum([chain.chain_mass.value_in_unit(unit.dalton)/c.Avogadro*unit.gram * n_copies 
                        for chain, n_copies in chain_info.items()], 0*unit.gram)
    
    target_volume = total_mass_g / target_density
    short_side_length = (target_volume / 6) ** (1/3)  # default long_side_scale_factor=6
    short_side_length = short_side_length.in_units_of(unit.nanometer) / unit.nanometer

    return short_side_length * np.array([[6, 0, 0], 
                                         [0, 1, 0], 
                                         [0, 0, 1]]) * unit.nanometer


# ====== Model building ======

def build_model(chain_info, target_box_vectors, use_random_placement=True, clash_cutoff=0.75, verbose=False):
    """
    Build a Modeller object with the given chain information and box vectors.
    
    This function places chains in the box using CALVADOS-style random placement:
    - If use_random_placement=True: use random placement with clash detection
    - If use_random_placement=False: center all chains at their original positions
    
    Args:
        chain_info (dict): 
            A dictionary with keys = chain objects and values = number of copies to place.
            Each chain object must have:
              - topology (OpenMM Topology)
              - min_rg_coords (initial coords)
        target_box_vectors (Quantity): 
            A 3x3 array of box vectors (in nm).
        use_random_placement (bool): 
            If True, use random placement with clash detection (CALVADOS style).
            If False, center each chain at its original position.
        clash_cutoff (float): 
            Minimum distance (nm) between atoms of different molecules. Default 0.75 nm.
        verbose (bool): 
            Print verbose output during placement.

    Returns:
        Modeller: 
            A Modeller object containing the combined system of all chains, 
            set with appropriately extended periodic box vectors.
    """
    print('Initializing the model...', flush=True)
    
    # Get box dimensions
    box_vectors = target_box_vectors.value_in_unit(unit.nanometer)
    Lx, Ly, Lz = box_vectors[0, 0], box_vectors[1, 1], box_vectors[2, 2]
    box = np.array([Lx, Ly, Lz])
    
    # Build chain list (expand copies)
    chain_list = []
    for chain, n_copies in chain_info.items():
        chain_list.extend([(chain, copy_idx) for copy_idx in range(n_copies)])
    
    total_chains = len(chain_list)
    print(f'  Total chains to place: {total_chains}', flush=True)
    
    if total_chains == 0:
        raise ValueError("No chains to place in chain_info")
    
    if use_random_placement:
        # ====== CALVADOS-style random placement ======
        print(f'  Using CALVADOS-style random placement (clash cutoff: {clash_cutoff} nm)', flush=True)
        
        # Start with first chain at center of box
        first_chain, first_copy_idx = chain_list[0]
        
        # Center the first molecule at the box center
        xinit = first_chain.min_rg_coords
        x0 = box / 2  # Box center
        first_positions = x0 + xinit
        
        model = app.Modeller(first_chain.topology, first_positions * unit.nanometer)
        
        # Track all positions for clash detection
        all_positions = [first_positions]
        
        # Place remaining chains
        for chain_idx in range(1, total_chains):
            chain_obj, copy_idx = chain_list[chain_idx]
            xinit = chain_obj.min_rg_coords
            
            if verbose:
                print(f'  Placing chain {chain_idx + 1}/{total_chains} ({chain_obj.chain_id} copy {copy_idx})...', flush=True)
            
            # Try random placement
            try:
                final_positions = random_placement(
                    box, all_positions, xinit, 
                    ntries=10000, cutoff=clash_cutoff, verbose=verbose
                )
            except ValueError as e:
                print(f'  WARNING: {e}', flush=True)
                print(f'  Falling back to center placement...', flush=True)
                final_positions = box / 2 + xinit
            
            model.add(chain_obj.topology, final_positions * unit.nanometer)
            all_positions.append(final_positions)
    
    else:
        # ====== Simple center placement (original behavior) ======
        print('  Using center placement (all chains at original positions)', flush=True)
        
        first_chain, first_copy_idx = chain_list[0]
        model = app.Modeller(first_chain.topology, first_chain.min_rg_coords)
        
        for chain_idx in range(1, total_chains):
            chain_obj, copy_idx = chain_list[chain_idx]
            model.add(chain_obj.topology, chain_obj.min_rg_coords * unit.nanometer)

    # Set periodic box vectors
    model.topology.setPeriodicBoxVectors(target_box_vectors)
    
    print('The model is built.', flush=True)
    print('The topology of the model:', model.topology)
    
    return model


# ====== System building ======

def build_mpipi_recharged_system_from_chains(chain_info,
                                              box_size=None,
                                              T=280*unit.kelvin,
                                              csx=150,
                                              CM_remover=True,
                                              periodic=True,
                                              use_random_placement=True,
                                              clash_cutoff=0.75,
                                              verbose=False):
    """
    Build an OpenMM System using Mpipi-Recharged forcefield from chain objects.
    
    This function:
    1. Calls get_compact_model() on each chain to relax and calculate Rg
    2. Builds a combined model using build_model() with CALVADOS-style random placement
    3. Adds forces using get_mpipi_system()
    
    Args:
        chain_info (dict):
            A dictionary with keys = chain objects (MDP, IDP, RNA, etc.) and values = number of copies.
            Each chain object must have:
              - topology (OpenMM Topology)
              - initial_coords (initial coordinates)
        box_size (list, optional):
            Box size [Lx, Ly, Lz] in nm. If provided, use this directly.
        T (Quantity, optional):
            Temperature in Kelvin (default 280 K).
        csx (float, optional):
            Ionic strength in mM (default 150).
        CM_remover (bool, optional):
            Whether to add center-of-mass motion remover (default True).
        periodic (bool, optional):
            Whether to use periodic boundary conditions (default True).
        use_random_placement (bool, optional):
            Whether to use CALVADOS-style random placement (default True).
        clash_cutoff (float, optional):
            Minimum distance (nm) between atoms of different molecules (default 0.75 nm).
        verbose (bool, optional):
            Print verbose output during placement.
    
    Returns:
        tuple: (system, model) where:
            - system: An OpenMM System object ready for simulation
            - model: A Modeller object containing the combined topology and positions
    """
    print('Relaxing monomers...', flush=True)
    
    # Step 1: Relax each chain and calculate Rg
    for chain in chain_info.keys():
        print(f'  Relaxing {chain.chain_id}...', flush=True)
        chain.get_compact_model(simulation_time=10*unit.nanosecond)
    
    # Step 2: Calculate target box vectors
    print('Calculating target box vectors...', flush=True)
    target_box_vectors = calculate_target_box_vectors(chain_info, box_size=box_size)
    print(f'  Box size: {target_box_vectors[0,0]} x {target_box_vectors[1,1]} x {target_box_vectors[2,2]}')
    
    # Step 3: Build combined model
    print('Building combined model...', flush=True)
    
    total_chains = sum(chain_info.values())
    if use_random_placement and total_chains > 1:
        model = build_model(
            chain_info, target_box_vectors, 
            use_random_placement=True, 
            clash_cutoff=clash_cutoff,
            verbose=verbose
        )
    else:
        print(f'  Placing {total_chains} chain(s) at center...', flush=True)
        model = build_model(
            chain_info, target_box_vectors, 
            use_random_placement=False
        )
    
    # Step 4: Build system with forces
    print('Building Mpipi-Recharged system...', flush=True)
    
    # Build globular_indices_dict
    globular_indices_dict = {
        chain.chain_id: chain.globular_indices for chain in chain_info.keys()
    }
    
    system = get_mpipi_system(
        np.array(model.positions),
        model.topology,
        globular_indices_dict,
        T.value_in_unit(unit.kelvin),
        csx,
        CM_remover=CM_remover,
        periodic=periodic
    )
    
    print('System built successfully.', flush=True)
    return system, model


# ====== Equilibration functions ======

def equilibrate_slab(model, 
                     target_box_vectors, 
                     chain_info, 
                     T=280*unit.kelvin, 
                     csx=150, 
                     pulling_time=20*unit.nanosecond, 
                     equi_time=400*unit.nanosecond):
    """
    Perform a two-stage equilibration:
    
    1) Pull all chains toward the center of the box to remove large voids (a "slab" pulling).
    2) Remove the pulling force, set the final box vectors, and equilibrate further.

    Args:
        model (Modeller): 
            A Modeller object with the combined system and initial positions.
        target_box_vectors (Quantity): 
            A 3x3 array of box vectors (in nm).
        chain_info (dict): 
            A dictionary with keys = chain objects and values = number of copies. 
            Used to retrieve `globular_indices`.
        T (Quantity, optional): 
            Temperature to run the simulation, default 280 K.
        csx (float, optional): 
            Ionic strength in mM, used to compute Debye length internally.
        pulling_time (Quantity, optional): 
            Duration of the pulling simulation (default 20 ns).
        equi_time (Quantity, optional): 
            Duration of the equilibration simulation after pulling (default 400 ns).
    """
    print('Setting up integrator...', flush=True)
    integrator = mm.LangevinMiddleIntegrator(T, 0.01/unit.picosecond, 10*unit.femtosecond)

    # Prepare the globular indices dictionary
    globular_indices_dict = {
        chain.chain_id: chain.globular_indices for chain in chain_info.keys()
    }

    print('Building mpipi system...', flush=True)
    # Calls your second-code function get_mpipi_system
    system = get_mpipi_system(
        np.array(model.positions), 
        model.topology, 
        globular_indices_dict, 
        T.value_in_unit(unit.kelvin), 
        csx, 
        CM_remover=False, 
        periodic=True
    )

    # Create a gentle pulling force to drive everything toward the box center
    pulling_force = mm.CustomExternalForce(
        'k*periodicdistance(x, y, z, x0, y0, z0)^2'
    )
    midpoint_x = 0.5*target_box_vectors[0][0]
    midpoint_y = 0.5*target_box_vectors[1][1]
    midpoint_z = 0.5*target_box_vectors[2][2]
    
    pulling_force.addGlobalParameter('x0', midpoint_x)
    pulling_force.addGlobalParameter('y0', midpoint_y)
    pulling_force.addGlobalParameter('z0', midpoint_z)
    pulling_force.addPerParticleParameter('k')

    # Add a small spring constant for each atom
    for atom in model.topology.atoms():
        pulling_force.addParticle(atom.index, [0.001])
    system.addForce(pulling_force)
    print('Pulling force is set.', flush=True)

    # Setting up the Simulation
    print('Initializing simulation...', flush=True)
    simulation = app.Simulation(
        model.topology, 
        system, 
        integrator, 
        platform=PLATFORM, 
        platformProperties=PROPERTIES
    )
    simulation.context.setPositions(model.positions)
    simulation.context.setPeriodicBoxVectors(*model.topology.getPeriodicBoxVectors())

    # Minimize
    simulation.minimizeEnergy()
    print('Energy minimized.', flush=True)

    # Pulling stage
    simulation.reporters.append(
        app.StateDataReporter(
            './output_pulling.dat', 10000, 
            step=True, potentialEnergy=True, 
            temperature=True, density=True, elapsedTime=True
        )
    )
    simulation.reporters.append(app.XTCReporter('./traj_pulling.xtc', 10000))
    print('Beginning pulling...', flush=True)
    simulation.step(int(pulling_time/(10*unit.femtosecond)))
    print('Pulling complete.', flush=True)

    # Remove pulling force
    system.removeForce(system.getNumForces()-1)
    print('Pulling force removed.', flush=True)

    # Update box vectors
    simulation.context.setPeriodicBoxVectors(*target_box_vectors)
    print('Box vectors updated.', flush=True)
    
    # Re-initialize the context so that box changes take effect
    state = simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    simulation.context.reinitialize()
    simulation.context.setState(state)
    print('Context reinitialized.', flush=True)

    # Minimize again
    model.topology.setPeriodicBoxVectors(target_box_vectors)
    simulation.minimizeEnergy()

    # Clear old reporters
    simulation.reporters = []

    # Equilibration stage
    simulation.reporters.append(
        app.StateDataReporter(
            './output_equi.dat', 10000, 
            step=True, potentialEnergy=True, 
            temperature=True, elapsedTime=True
        )
    )
    simulation.reporters.append(app.XTCReporter('./traj_equi.xtc', 50000))
    print('Beginning equilibration...', flush=True)
    simulation.step(int(equi_time/(10*unit.femtosecond)))
    
    # Final state
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    app.PDBFile.writeFile(model.topology, positions, open('./equi_model.pdb', 'w'))
    simulation.saveState('equi_state.xml')
    print('Equilibration complete.', flush=True)


def build_and_equilibrate_model(chain_info, 
                                box_size=None,
                                target_density=0.1*unit.gram/unit.centimeter**3, 
                                T=280*unit.kelvin, 
                                csx=150,
                                pulling_time=20*unit.nanosecond, 
                                equi_time=1000*unit.nanosecond,
                                use_random_placement=True,
                                clash_cutoff=0.75,
                                verbose=False):
    """
    Wrapper function to:
      1) Relax each chain monomer in isolation.
      2) Calculate and build the combined model with the desired box size (CALVADOS-style).
      3) Perform slab-pulling followed by a full equilibration.

    Args:
        chain_info (dict):
            A dictionary with keys = chain objects and values = number of copies. 
            Each chain object must have:
              - get_compact_model(simulation_time): a method to relax the chain individually
              - chain_mass, min_rg_coords, max_rg, min_rg
        box_size (list, optional):
            Box size [Lx, Ly, Lz] in nm. If provided, use this directly.
        target_density (Quantity, optional):
            Target density in g/cm^3.
        T (Quantity, optional):
            Simulation temperature.
        csx (float, optional):
            Ionic strength in mM (for Debye length calculations).
        pulling_time (Quantity, optional):
            Duration of the slab-pulling simulation.
        equi_time (Quantity, optional):
            Duration of the final equilibration simulation.
        use_random_placement (bool, optional):
            Whether to use CALVADOS-style random placement (default True).
        clash_cutoff (float, optional):
            Minimum distance (nm) between atoms of different molecules (default 0.75 nm).
        verbose (bool, optional):
            Print verbose output during placement.
    """
    
    print('Relaxing monomers...')
    for chain in chain_info.keys():
        chain.get_compact_model(simulation_time=10*unit.nanosecond)
    
    # Build combined model
    target_box_vectors = calculate_target_box_vectors(chain_info, 
                                                      box_size=box_size,
                                                      target_density=target_density)
    
    model = build_model(
        chain_info, target_box_vectors,
        use_random_placement=use_random_placement,
        clash_cutoff=clash_cutoff,
        verbose=verbose
    )

    # Equilibrate
    equilibrate_slab(model, target_box_vectors, chain_info, 
                     T=T, csx=csx, 
                     pulling_time=pulling_time, 
                     equi_time=equi_time)
