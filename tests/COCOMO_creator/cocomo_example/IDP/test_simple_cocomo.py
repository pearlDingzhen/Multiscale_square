#!/usr/bin/env python3
"""
Simplified test to verify COCOMO system building.
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


def test_system_building():
    """Test that new COCOMO can build a valid system."""
    pdb_file = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/IDP/4chain.pdb'
    box_size = 30.0
    
    # Read PDB with OpenMM to get chain info
    pdb = PDBFile(pdb_file)
    
    # Count chains from topology
    n_chains = len(list(pdb.topology.chains()))
    print(f"Number of chains in PDB: {n_chains}")
    
    # Count residues
    n_residues = pdb.topology._numResidues
    print(f"Number of residues: {n_residues}")
    
    # Count atoms
    n_atoms = pdb.topology._numAtoms
    print(f"Number of atoms: {n_atoms}")
    
    # Build chain_ids list using topology atoms
    atoms_list = list(pdb.topology.atoms())
    residues_list = list(pdb.topology.residues())
    chains_list = list(pdb.topology.chains())
    
    # Build chain_id for each residue
    # Map atom index to chain id
    atom_chain_map = {}
    for chain in chains_list:
        for atom in chain.atoms():
            atom_chain_map[atom.index] = chain.id
    
    chain_ids = []
    for res in residues_list:
        # Get first atom in residue to determine chain
        res_atoms = list(res.atoms())
        if res_atoms:
            chain_ids.append(atom_chain_map.get(res_atoms[0].index, 'A'))
        else:
            chain_ids.append('A')
    
    print(f"Chain IDs: {set(chain_ids)}")
    print(f"Chain ID list length: {len(chain_ids)}")
    
    # Build global sequence from topology
    aa_3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    residues_list = list(pdb.topology.residues())
    global_sequence = ''.join(aa_3to1.get(res.name, 'A') for res in residues_list)
    print(f"Sequence length: {len(global_sequence)}")
    print(f"Sequence (first 50): {global_sequence[:50]}")
    
    # Get positions
    positions = np.array(pdb.positions.value_in_unit(nanometers))
    print(f"Positions shape: {positions.shape}")
    
    # Build topology_info
    # Map chain IDs (letters) to integers
    unique_chains = sorted(set(chain_ids))
    chain_id_to_int = {c: i+1 for i, c in enumerate(unique_chains)}
    int_chain_ids = [chain_id_to_int[c] for c in chain_ids]
    
    topology_info = {
        'global_sequence': global_sequence,
        'chain_ids': int_chain_ids,
        'folded_domains': [0] * len(global_sequence),
        'component_names': ['chain'] * len(global_sequence),
        'local_residue_indices': list(range(1, len(global_sequence) + 1))
    }
    
    # Create COCOMO system
    print("\nCreating COCOMO system...")
    cocomo = COCOMO(
        box_size=box_size,
        topology_info=topology_info,
        positions=positions,
        surf=0.7,
        resources='CPU'
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
    
    print("\n✓ SUCCESS: COCOMO system built and energy computed successfully!")
    

if __name__ == '__main__':
    test_system_building()

