#!/usr/bin/env python
"""
Example usage of ms2_cg2all for converting CG structures to all-atom models.

This demonstrates how to use the integrated cg2all functionality within multiscale2.
"""

import sys
sys.path.insert(0, '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square')

from multiscale2.extern.ms2_cg2all import convert_cg2all

# Example 1: Convert a Calpha-based model to all-atom structure
print("Example 1: Calpha to all-atom conversion")
convert_cg2all(
    in_pdb_fn="input_structure.ca.pdb",
    out_fn="output_structure.all.pdb",
    model_type="CalphaBasedModel"
)
print("✓ Conversion complete!")

# Example 2: Convert with GPU acceleration (if available)
print("\nExample 2: Using GPU acceleration")
convert_cg2all(
    in_pdb_fn="input_structure.ca.pdb",
    out_fn="output_structure.all.pdb",
    model_type="CalphaBasedModel",
    device="cuda"  # or "cpu" for CPU-only execution
)
print("✓ GPU conversion complete!")

# Example 3: Convert with atom preservation
print("\nExample 3: Preserving CA coordinates from input")
convert_cg2all(
    in_pdb_fn="input_structure.ca.pdb",
    out_fn="output_structure.all.pdb",
    model_type="CalphaBasedModel",
    fix_atom=True  # Keep CA positions from input
)
print("✓ Conversion with atom preservation complete!")

print("\n" + "="*60)
print("All examples completed successfully!")
print("="*60)





