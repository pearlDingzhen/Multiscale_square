#!/usr/bin/env python3
"""
Test script to compare old and new COCOMO implementations for MDP (folded proteins).
Uses a PDB with proper segment_id for OLD implementation compatibility.
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

from cocomo2_creator import COCOMO


def run_old_cocomo(pdb_file, surface_file, box_size, domain, count):
    """
    Run OLD COCOMO implementation as subprocess and parse energy output.
    
    Returns:
        dict: energy breakdown by force name
    """
    # Create temp directory with necessary files
    import shutil
    import tempfile
    temp_dir = tempfile.mkdtemp()
    
    # Copy PDB and surface files
    shutil.copy(pdb_file, temp_dir + '/system_with_segid.pdb')
    shutil.copy(surface_file, temp_dir + '/surface')
    
    # Build command matching: python cocomo2_old.py -n 16 -d 3 76 106 176 192 260 320 334 -b 30 -f surface -p . -f system_with_segid.pdb
    cmd = [
        sys.executable, 'cocomo2_old.py',
        '-p', temp_dir + '/',           # Path to working directory
        '-n', str(count),               # Number of chains/monomers
        '-b', str(box_size),            # Box size in nm
        '-f', 'surface',                # Surface file name
        '-f', 'system_with_segid.pdb'   # PDB file name (appears twice in command)
    ]
    
    # Add domain parameters (-d flag with all boundaries in one flag)
    if domain and len(domain) >= 2:
        cmd.append('-d')
        for d in domain:
            cmd.append(str(d))
    
    # Run from the MDP directory where cocomo2_old.py is located
    mdp_dir = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP'
    result = subprocess.run(
        cmd,
        cwd=mdp_dir,
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


def build_system_new(pdb_file, box_size, surface_file, domain, surf=0.7):
    """Build system using new COCOMO class."""
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
    
    # Keep positions as Quantity (don't convert to plain numpy array)
    positions_array = pdb.positions.value_in_unit(nanometers)
    positions = np.array([[float(x), float(y), float(z)] * nanometers 
                          for x, y, z in positions_array])
    
    # Load surface values from file
    surface_values = np.loadtxt(surface_file)
    n_chains = len(unique_chains)
    
    # CRITICAL: Match OLD COCOMO's surface_vector processing
    # OLD COCOMO sets NON-domain regions to 5.0, keeps domain regions as original surface values
    # This affects the SASA scaling factor calculation!
    if domain and len(domain) >= 2:
        # Make a copy to modify
        surface_values_processed = surface_values.copy()
        
        # Process domain boundaries (1-based indexing in OLD COCOMO)
        # Set non-domain regions to 5.0, keep domain regions as original
        for j in range(len(domain)):
            if j == 0 and domain[j] > 1:
                # Before first domain: set to 5.0
                surface_values_processed[:domain[j] - 1] = 5.0
            elif j % 2 == 0:
                # Between domains: set to 5.0
                surface_values_processed[domain[j - 1] - 1:domain[j] - 1] = 5.0
            elif j == len(domain) - 1 and domain[j] < len(surface_values_processed):
                # After last domain: set to 5.0
                surface_values_processed[domain[j] - 1:] = 5.0
        
        print(f"  Processed surface values: domain regions keep original values, non-domain set to 5.0")
        print(f"  Domain regions: {[(domain[k], domain[k+1]) for k in range(0, len(domain), 2)]}")
        
        # Build folded_domains for elastic network (1.0 = folded, 5.0 = unfolded)
        folded_single = np.ones(len(surface_values)) * 5.0
        for k in range(0, len(domain), 2):
            start_idx = domain[k] - 1
            end_idx = domain[k + 1]
            if start_idx < len(folded_single) and end_idx <= len(folded_single):
                folded_single[start_idx:end_idx] = 1.0
    else:
        # No domains specified - treat all as unfolded (5.0) for surface, folded for ENM
        surface_values_processed = np.ones(len(surface_values)) * 5.0
        folded_single = np.ones(len(surface_values)) * 5.0
        print("  No domains specified, treating all as unfolded (5.0)")
    
    # Tile for all chains (match OLD COCOMO: surface_vector.tolist() * count)
    folded_domains = np.tile(folded_single, n_chains)
    sasa_values = np.tile(surface_values_processed, n_chains)
    
    topology_info = {
        'global_sequence': global_sequence,
        'chain_ids': int_chain_ids,
        'folded_domains': folded_domains.tolist(),
        'component_names': ['TDP43'] * len(global_sequence),
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


def get_force_energies(system, positions, platform_name='CPU', debug_sasa=False):
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
    # Force order in OLD COCOMO: Bond, Angle, CMMotionRemover, Electrostatic, VdW, Cation-pi, Pi-pi, ElasticNetwork
    # Force order in NEW COCOMO: Bond, Angle, CMMotionRemover, Electrostatic, VdW, Cation-pi, Pi-pi, ElasticNetwork
    # Both implementations have the same force order
    force_names = {
        0: 'HarmonicBondForce',
        1: 'HarmonicAngleForce',
        2: 'CMMotionRemover',
        3: 'Electrostatic',
        4: 'Van der Waals',
        5: 'Cation-pi',
        6: 'Pi-pi',
        7: 'ElasticNetwork'
    }
    
    # Debug: Check SASA scaling factors for VdW force
    if debug_sasa and system.getNumForces() > 4:
        vdw_force = system.getForce(4)  # Van der Waals force
        if hasattr(vdw_force, 'getNumParticles'):
            print(f"\n  Debug: VdW force has {vdw_force.getNumParticles()} particles")
            # Sample first 10 particles' S parameters
            for i in range(min(10, vdw_force.getNumParticles())):
                params = vdw_force.getParticleParameters(i)
                if len(params) >= 3:
                    print(f"    Particle {i}: S={params[2]:.6f}")
    
    for i in range(system.getNumForces()):
        state = context.getState(getEnergy=True, groups={i})
        energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        name = force_names.get(i, f'Force {i}')
        force_energies[name] = energy
    
    force_energies['Particles'] = system.getNumParticles()
    force_energies['Forces'] = system.getNumForces()
    
    return force_energies


def main():
    """Main test function for MDP comparison."""
    # Use the new PDB with segment_id (same file for both OLD and NEW COCOMO)
    pdb_file_with_segid = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/system_with_segid.pdb'
    surface_file = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/surface'
    
    # Match the command: python cocomo2_old.py -n 16 -d 3 76 106 176 192 260 320 334 -b 30 -f surface -p . -f system_with_segid.pdb
    box_size = 40.0  # -b 30
    count = 16       # -n 16 (number of monomers)
    
    # Domain boundaries: -d 3 76 106 176 192 260 320 334
    # These define the folded regions: [3,76], [106,176], [192,260], [320,334]
    domain = [3, 76, 106, 176, 192, 260, 320, 334]
    
    # For OLD COCOMO: 
    # - It expects a single-chain PDB and replicates it count times
    # - But system_with_segid.pdb already has multiple chains with segment_id
    # - OLD's chain detection uses segment_id to determine chain boundaries
    # - So we use the same PDB file and let OLD handle the topology
    
    # For NEW COCOMO:
    # - We need to extract chain information from the PDB's chain IDs
    # - Or use the system_with_segid.pdb which has proper chain structure
    
    print("=" * 70)
    print("Comparing Old vs New COCOMO Implementation for MDP")
    print("=" * 70)
    print(f"PDB file: {pdb_file_with_segid}")
    print(f"Surface: {surface_file}")
    print(f"Box: {box_size} nm")
    print(f"Monomer count: {count}")
    print(f"Domains: {domain}")
    print(f"  - Domain 1: [3, 76]")
    print(f"  - Domain 2: [106, 176]")
    print(f"  - Domain 3: [192, 260]")
    print(f"  - Domain 4: [320, 334]")
    
    # Run OLD COCOMO
    print("\n" + "=" * 70)
    print("Running OLD COCOMO...")
    print("=" * 70)
    old_energies = run_old_cocomo(pdb_file_with_segid, surface_file, box_size, domain, count)
    
    # Build NEW COCOMO using the same PDB file
    print("\n" + "=" * 70)
    print("Building NEW COCOMO...")
    print("=" * 70)
    system_new, top_new = build_system_new(pdb_file_with_segid, box_size, surface_file, domain, surf=0.7)
    
    # Get positions with units for NEW COCOMO energy calculation
    # Use the same PDB file as OLD COCOMO to ensure identical coordinates
    pdb_for_new = PDBFile(pdb_file_with_segid)
    positions_array = pdb_for_new.positions.value_in_unit(nanometers)
    positions_with_units = np.array([[float(x), float(y), float(z)] * nanometers 
                                     for x, y, z in positions_array])
    print(f"Using positions from: {pdb_file_with_segid}")
    print(f"Number of positions: {len(positions_with_units)}")
    new_energies = get_force_energies(system_new, positions_with_units, debug_sasa=True)
    
    # Compare
    print("\n" + "=" * 70)
    print("Energy Comparison")
    print("=" * 70)
    print(f"{'Force':<25} | {'Old (kJ/mol)':<15} | {'New (kJ/mol)':<15} | {'Ratio':<10}")
    print("-" * 70)
    
    # Compare common forces
    common_forces = ['CMMotionRemover', 'HarmonicBondForce', 'HarmonicAngleForce', 
                     'Electrostatic', 'Van der Waals', 'Cation-pi', 'Pi-pi', 'ElasticNetwork']
    
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
