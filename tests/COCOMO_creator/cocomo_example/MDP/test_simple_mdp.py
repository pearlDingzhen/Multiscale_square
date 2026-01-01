#!/usr/bin/env python3
"""
Simplified test to verify COCOMO system building for MDP (folded proteins).
Uses system.pdb which contains 16 TDP43 molecules (6624 residues total).
"""
import numpy as np
from openmm.unit import *
from openmm import *
from openmm.app import *
import os
import sys

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
sys.path.insert(0, project_root)

from cocomo2_creator import COCOMO


def test_mdp_system_building():
    """Test that new COCOMO can build a valid MDP system."""
    pdb_file = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/system.pdb'
    box_size = 40.0  # nm (from config_mdp.yaml)
    
    # Read PDB with OpenMM
    pdb = PDBFile(pdb_file)
    
    # Count chains, residues, atoms from topology
    n_chains = len(list(pdb.topology.chains()))
    n_residues = pdb.topology._numResidues
    n_atoms = pdb.topology._numAtoms
    
    print(f"Number of chains in PDB: {n_chains}")
    print(f"Number of residues: {n_residues}")
    print(f"Number of atoms: {n_atoms}")
    
    # Build chain_id for each residue using topology atoms
    atoms_list = list(pdb.topology.atoms())
    residues_list = list(pdb.topology.residues())
    chains_list = list(pdb.topology.chains())
    
    # Map atom index to chain id
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
    
    print(f"Unique chains: {sorted(set(chain_ids))}")
    print(f"Chain ID list length: {len(chain_ids)}")
    
    # Build global sequence from topology
    aa_3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    global_sequence = ''.join(aa_3to1.get(res.name, 'A') for res in residues_list)
    print(f"Sequence length: {len(global_sequence)}")
    print(f"Sequence (first 50): {global_sequence[:50]}")
    print(f"Sequence (last 50): {global_sequence[-50:]}")
    
    # Get positions
    positions = np.array(pdb.positions.value_in_unit(nanometers))
    print(f"Positions shape: {positions.shape}")
    
    # Map chain IDs to integers (1-based)
    unique_chains = sorted(set(chain_ids))
    chain_id_to_int = {c: i+1 for i, c in enumerate(unique_chains)}
    int_chain_ids = [chain_id_to_int[c] for c in chain_ids]
    
    # Load SASA values from surface file
    # This is what OLD implementation uses
    surface_values = np.loadtxt('/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/surface')
    
    print(f"\nSASA values:")
    print(f"  - Source: surface file (414 values)")
    print(f"  - Non-5.0 (folded): {np.sum(surface_values != 5.0)}")
    print(f"  - == 5.0 (unfolded): {np.sum(surface_values == 5.0)}")
    
    # Build topology_info
    topology_info = {
        'global_sequence': global_sequence,
        'chain_ids': int_chain_ids,
        'folded_domains': [1 if v != 5.0 else 0 for v in surface_values] * 16,  # Repeat for 16 molecules
        'component_names': ['TDP43'] * len(global_sequence),
        'local_residue_indices': list(range(1, len(global_sequence) + 1)),
        'sasa_values': np.tile(surface_values, 16)  # Repeat for 16 molecules
    }
    
    print(f"\nTopology info:")
    print(f"  - global_sequence: {len(topology_info['global_sequence'])} residues")
    print(f"  - chain_ids: {len(topology_info['chain_ids'])} values")
    print(f"  - folded_domains: {sum(topology_info['folded_domains'])} folded / {len(topology_info['folded_domains']) - sum(topology_info['folded_domains'])} unfolded")
    print(f"  - sasa_values: {len(topology_info['sasa_values'])} values")
    
    # Create COCOMO system
    print("\nCreating COCOMO system...")
    cocomo = COCOMO(
        box_size=box_size,
        topology_info=topology_info,
        positions=positions,
        surf=0.7,
        resources='CPU',
        sasa_values=topology_info['sasa_values']
    )
    
    print("Building system...")
    system = cocomo.create_system()
    
    print(f"\n=== System Info ===")
    print(f"Particles: {system.getNumParticles()}")
    print(f"Forces: {system.getNumForces()}")
    
    # List forces
    for i in range(system.getNumForces()):
        print(f"  Force {i}: {type(system.getForce(i)).__name__}")
    
    # Test energy calculation
    print("\nComputing energy...")
    topology = cocomo._build_topology()
    
    integrator = LangevinIntegrator(298 * kelvin, 0.01 / picosecond, 0.01 * picoseconds)
    platform = Platform.getPlatformByName('CPU')
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Initial energy: {energy}")
    
    print("\n✓ SUCCESS: MDP COCOMO system built and energy computed successfully!")
    
    # Additional info
    print(f"\n=== Additional Info ===")
    print(f"Box size: {box_size} nm")
    print(f"Number of TDP43 molecules: {n_chains}")
    print(f"Residues per TDP43: {n_residues // n_chains}")
    print(f"Surface scaling (surf): 0.7")


if __name__ == '__main__':
    test_mdp_system_building()




