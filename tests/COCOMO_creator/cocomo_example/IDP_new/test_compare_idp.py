#!/usr/bin/env python3
"""
Test script to compare old and new COCOMO implementations for IDP (intrinsically disordered proteins).
Uses a PDB with proper segment_id for OLD implementation compatibility.
IDP systems have no domain boundaries, so all surface values are set to 5.0.
"""
import subprocess
import sys
import re
import numpy as np
from openmm.unit import *
from openmm import *
from openmm.app import *
import os

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
sys.path.insert(0, project_root)

# Import COCOMO from MDP directory (same implementation for both MDP and IDP)
mdp_dir = os.path.dirname(os.path.dirname(__file__)) + '/MDP'
sys.path.insert(0, mdp_dir)
from cocomo2_creator import COCOMO


def run_old_cocomo(pdb_file, box_size, count):
    """
    Run OLD COCOMO implementation as subprocess and parse energy output.
    IDP systems: no domain parameters, no surface file.
    
    Returns:
        dict: energy breakdown by force name
    """
    # Create temp directory with necessary files
    import shutil
    import tempfile
    temp_dir = tempfile.mkdtemp()
    
    # Copy PDB file
    shutil.copy(pdb_file, temp_dir + '/4chain_with_segid.pdb')
    
    # Build command matching: python cocomo2_old.py -n 4 -b 30 -f 4chain_with_segid.pdb -p .
    cmd = [
        sys.executable, 'cocomo2_old.py',
        '-p', temp_dir + '/',           # Path to working directory
        '-n', str(count),               # Number of chains/monomers
        '-b', str(box_size),            # Box size in nm
        '-f', '4chain_with_segid.pdb'   # PDB file name
    ]
    
    # IDP: No domain parameters (-d flag)
    # IDP: No surface file (-s flag), all surface values will be 5.0
    
    # Run from the IDP_new directory where cocomo2_old.py is located
    idp_dir = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/IDP_new'
    result = subprocess.run(
        cmd,
        cwd=idp_dir,
        capture_output=True,
        text=True
    )
    
    print(f"OLD COCOMO command: {' '.join(cmd)}")
    print(f"OLD COCOMO output:\n{result.stdout}")
    if result.returncode != 0:
        print(f"OLD COCOMO error:\n{result.stderr}")
    
    # Parse energy output using robust regex
    energies = {}
    for line in result.stdout.split('\n'):
        # Match pattern like 'Force Name | Energy'
        match = re.match(r'^(.+?)\s*\|\s*([\d.-]+)', line.strip())
        if match:
            name = match.group(1).strip()
            energy = float(match.group(2))
            # Normalize force names to match NEW COCOMO naming
            name_mapping = {
                'Electrostatic (Nonbonded)': 'Electrostatic',
                'Van der Waals (LJ)': 'Van der Waals',
            }
            name = name_mapping.get(name, name)
            energies[name] = energy
    
    # Also get totals
    for line in result.stdout.split('\n'):
        if 'Particles:' in line:
            match = re.search(r'Particles:\s*(\d+)', line)
            if match:
                energies['Particles'] = int(match.group(1))
        if 'Forces:' in line:
            match = re.search(r'Forces:\s*(\d+)', line)
            if match:
                energies['Forces'] = int(match.group(1))
    
    return energies


def build_system_new(pdb_file, box_size, surf=0.7):
    """Build system using new COCOMO class for IDP."""
    from biopandas.pdb import PandasPdb
    
    pdb = PDBFile(pdb_file)
    
    # Build chain_ids list from OpenMM topology chains
    residues_list = list(pdb.topology.residues())
    chains_list = list(pdb.topology.chains())
    
    atom_chain_map = {}
    for chain in chains_list:
        for atom in chain.atoms():
            atom_chain_map[atom.index] = chain.id
    
    chain_ids = []
    for res in residues_list:
        res_atoms = list(res.atoms())
        if res_atoms:
            chain_ids.append(atom_chain_map.get(res_atoms[0].index, 'A'))
        else:
            chain_ids.append('A')
    
    # Map chain IDs to integers (1-based)
    unique_chains = sorted(set(chain_ids))
    chain_id_to_int = {c: i+1 for i, c in enumerate(unique_chains)}
    int_chain_ids = [chain_id_to_int[c] for c in chain_ids]
    
    # Build global sequence
    aa_3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    global_sequence = ''.join(aa_3to1.get(res.name, 'A') for res in residues_list)
    
    # Keep positions as Quantity
    positions_array = pdb.positions.value_in_unit(nanometers)
    positions = np.array([[float(x), float(y), float(z)] * nanometers 
                          for x, y, z in positions_array])
    
    # IDP: All surface values are 5.0 (unfolded)
    # OLD COCOMO: surface_vector = [5] * len(atom_list) * count
    n_chains = len(unique_chains)
    n_residues_per_chain = len(global_sequence) // n_chains
    
    # All residues are unfolded (5.0) for IDP
    surface_values_processed = np.ones(n_residues_per_chain) * 5.0
    folded_single = np.ones(n_residues_per_chain) * 5.0  # All unfolded
    
    print(f"  IDP system: All surface values set to 5.0 (unfolded)")
    print(f"  Residues per chain: {n_residues_per_chain}")
    print(f"  Number of chains: {n_chains}")
    
    # Tile for all chains (match OLD COCOMO: [5] * len(atom_list) * count)
    folded_domains = np.tile(folded_single, n_chains)
    sasa_values = np.tile(surface_values_processed, n_chains)
    
    topology_info = {
        'global_sequence': global_sequence,
        'chain_ids': int_chain_ids,
        'folded_domains': folded_domains.tolist(),
        'component_names': ['TDP43_CTD'] * len(global_sequence),
        'local_residue_indices': list(range(1, len(global_sequence) + 1)),
        'sasa_values': sasa_values
    }
    
    cocomo = COCOMO(
        box_size=box_size,
        topology_info=topology_info,
        positions=positions,
        surf=surf,
        resources='CPU'
    )
    
    system = cocomo.create_system()
    return system, cocomo._build_topology()


def get_force_energies(system, positions, platform_name='CPU'):
    """Get each force's energy value."""
    for i, force in enumerate(system.getForces()):
        if i > 31:
            raise ValueError("OpenMM supports max 32 Force Groups")
        force.setForceGroup(i)
    
    integrator = VerletIntegrator(0.001)
    platform = Platform.getPlatformByName(platform_name)
    context = Context(system, integrator, platform)
    context.setPositions(positions)
    
    force_energies = {}
    # Force order in OLD COCOMO: Bond, Angle, CMMotionRemover, Electrostatic, VdW, Cation-pi, Pi-pi
    # IDP has no ElasticNetwork
    force_names = {
        0: 'HarmonicBondForce',
        1: 'HarmonicAngleForce',
        2: 'CMMotionRemover',
        3: 'Electrostatic',
        4: 'Van der Waals',
        5: 'Cation-pi',
        6: 'Pi-pi'
        # No ElasticNetwork for IDP
    }
    
    for i in range(system.getNumForces()):
        state = context.getState(getEnergy=True, groups={i})
        energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        name = force_names.get(i, f'Force {i}')
        force_energies[name] = energy
    
    force_energies['Particles'] = system.getNumParticles()
    force_energies['Forces'] = system.getNumForces()
    
    return force_energies


def main():
    """Main test function for IDP comparison."""
    # Use the PDB with segment_id (same file for both OLD and NEW COCOMO)
    pdb_file_with_segid = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/IDP_new/4chain_with_segid.pdb'
    
    # Match the command: python cocomo2_old.py -n 4 -b 30 -f 4chain_with_segid.pdb -p .
    box_size = 30.0  # -b 30
    count = 4        # -n 4 (number of chains)
    
    # IDP: No domain boundaries, no surface file
    # All surface values are 5.0 (unfolded)
    
    print("=" * 70)
    print("Comparing Old vs New COCOMO Implementation for IDP")
    print("=" * 70)
    print(f"PDB file: {pdb_file_with_segid}")
    print(f"Box: {box_size} nm")
    print(f"Monomer count: {count}")
    print(f"System type: IDP (no domains, all surface values = 5.0)")
    
    # Run OLD COCOMO
    print("\n" + "=" * 70)
    print("Running OLD COCOMO...")
    print("=" * 70)
    old_energies = run_old_cocomo(pdb_file_with_segid, box_size, count)
    
    # Build NEW COCOMO using the same PDB file
    print("\n" + "=" * 70)
    print("Building NEW COCOMO...")
    print("=" * 70)
    system_new, top_new = build_system_new(pdb_file_with_segid, box_size, surf=0.7)
    
    # Get positions with units for NEW COCOMO energy calculation
    # Use the same PDB file as OLD COCOMO to ensure identical coordinates
    pdb_for_new = PDBFile(pdb_file_with_segid)
    positions_array = pdb_for_new.positions.value_in_unit(nanometers)
    positions_with_units = np.array([[float(x), float(y), float(z)] * nanometers 
                                     for x, y, z in positions_array])
    print(f"Using positions from: {pdb_file_with_segid}")
    print(f"Number of positions: {len(positions_with_units)}")
    new_energies = get_force_energies(system_new, positions_with_units)
    
    # Compare
    print("\n" + "=" * 70)
    print("Energy Comparison")
    print("=" * 70)
    print(f"{'Force':<25} | {'Old (kJ/mol)':<15} | {'New (kJ/mol)':<15} | {'Ratio':<10}")
    print("-" * 70)
    
    # Compare common forces (IDP has no ElasticNetwork)
    common_forces = ['CMMotionRemover', 'HarmonicBondForce', 'HarmonicAngleForce', 
                     'Electrostatic', 'Van der Waals', 'Cation-pi', 'Pi-pi']
    
    all_match = True
    for name in common_forces:
        if name in old_energies and name in new_energies:
            e_old = old_energies[name]
            e_new = new_energies[name]
            
            if abs(e_old) > 0.001:
                ratio = abs(e_new / e_old) if e_old != 0 else float('inf')
            else:
                ratio = float('inf') if abs(e_new) > 0.001 else 1.0
            
            match_str = "✓" if ratio < 1.1 else "~" if ratio < 10 else "✗"
            if ratio >= 1.1:
                all_match = False
            
            print(f"{name:<25} | {e_old:>15.4f} | {e_new:>15.4f} | {ratio:>8.2f}x {match_str}")
    
    # Particles and Forces
    print(f"\nParticles: OLD={old_energies.get('Particles', 'N/A')}, NEW={new_energies.get('Particles', 'N/A')}")
    print(f"Forces: OLD={old_energies.get('Forces', 'N/A')}, NEW={new_energies.get('Forces', 'N/A')}")
    
    print("\n" + "=" * 70)
    if all_match:
        print("SUCCESS! Old and New COCOMO implementations match closely!")
    else:
        print("Some differences detected - see details above")
    print("=" * 70)


if __name__ == '__main__':
    main()

