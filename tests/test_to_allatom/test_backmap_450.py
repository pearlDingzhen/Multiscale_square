#!/usr/bin/env python
"""
Test script to backmap 450.pdb using CalphaBasedModel (without fix).
"""

import time
from multiscale2.extern.ms2_cg2all import convert_cg2all

# Input and output paths
input_pdb = "/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/test_to_allatom/450.pdb"
output_pdb = "/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/test_to_allatom/450.all.pdb"

print("Backmapping 450.pdb using CalphaBasedModel (without fix)...")
print(f"Input:  {input_pdb}")
print(f"Output: {output_pdb}")
print()

# Time the backmapping process
start_time = time.time()

convert_cg2all(
    in_pdb_fn=input_pdb,
    out_fn=output_pdb,
    model_type="CalphaBasedModel",
    fix_atom=False,
    device="cpu"  # Change to "cpu" if CUDA is not available
)

end_time = time.time()
elapsed_time = end_time - start_time

print()
print(f"Done!")
print(f"Backmapping completed in {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
