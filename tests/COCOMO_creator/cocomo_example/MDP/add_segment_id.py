#!/usr/bin/env python3
"""
Add segment_id to PDB file based on chain_id.
This allows OLD COCOMO implementation to properly identify chains.
"""
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np


def add_segment_id_to_pdb(input_pdb, output_pdb):
    """
    Read a PDB file and write a new one with segment_id set to chain_id.
    
    Args:
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path
    """
    universe = mda.Universe(input_pdb)
    
    # Create a new Universe with segment_id = chain_id
    # We'll manually write the PDB with proper segment_id
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    # Parse and modify ATOM records
    new_lines = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Keep first 72 characters (standard PDB columns)
            # segment_id is columns 73-76 (0-indexed: 72-75)
            chain_id = line[21]  # Chain ID is column 22 (0-indexed: 21)
            new_line = line[:72] + chain_id.ljust(4) + line[76:]
            new_lines.append(new_line)
        else:
            new_lines.append(line)
    
    with open(output_pdb, 'w') as f:
        f.writelines(new_lines)
    
    print(f"Created {output_pdb} with segment_id from chain_id")
    
    # Verify
    from biopandas.pdb import PandasPdb
    pdb = PandasPdb().read_pdb(output_pdb)
    print(f"segment_id unique values: {pdb.df['ATOM']['segment_id'].unique()}")


if __name__ == '__main__':
    input_file = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/system.pdb'
    output_file = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/system_with_segid.pdb'
    
    add_segment_id_to_pdb(input_file, output_file)




