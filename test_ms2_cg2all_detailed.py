#!/usr/bin/env python
"""
Detailed test script for ms2_cg2all integration.
"""

import sys
import os

# Add the project root to the path
sys.path.insert(0, '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square')

print("="*50)
print("Testing ms2_cg2all Integration")
print("="*50)

# Test 1: Import the lib module directly
print("\n1. Testing direct lib module imports...")
try:
    from multiscale2.extern.ms2_cg2all.lib import libconfig
    print("✓ Successfully imported libconfig")
except ImportError as e:
    print(f"✗ Failed to import libconfig: {e}")

# Test 2: Import residue_constants_base
print("\n2. Testing residue_constants_base...")
try:
    from multiscale2.extern.ms2_cg2all.lib.residue_constants_base import AMINO_ACID_s
    print(f"✓ Successfully imported residue_constants_base, AMINO_ACID_s has {len(AMINO_ACID_s)} items")
except ImportError as e:
    print(f"✗ Failed to import residue_constants_base: {e}")

# Test 3: Import residue_constants
print("\n3. Testing residue_constants...")
try:
    from multiscale2.extern.ms2_cg2all.lib.residue_constants import read_coarse_grained_topology
    print("✓ Successfully imported residue_constants")
except ImportError as e:
    print(f"✗ Failed to import residue_constants: {e}")

# Test 4: Import libpdb
print("\n4. Testing libpdb...")
try:
    from multiscale2.extern.ms2_cg2all.lib.libpdb import PDB
    print("✓ Successfully imported libpdb")
except ImportError as e:
    print(f"✗ Failed to import libpdb: {e}")

# Test 5: Import numpy_basics
print("\n5. Testing numpy_basics...")
try:
    from multiscale2.extern.ms2_cg2all.lib.numpy_basics import torsion_angle
    print("✓ Successfully imported numpy_basics")
except ImportError as e:
    print(f"✗ Failed to import numpy_basics: {e}")

# Test 6: Import torch_basics
print("\n6. Testing torch_basics...")
try:
    from multiscale2.extern.ms2_cg2all.lib.torch_basics import v_size
    print("✓ Successfully imported torch_basics")
except ImportError as e:
    print(f"✗ Failed to import torch_basics: {e}")

# Test 7: Import SE3Transformer
print("\n7. Testing SE3Transformer...")
try:
    from multiscale2.extern.ms2_se3transformer import Fiber, SE3Transformer
    print("✓ Successfully imported SE3Transformer")
except ImportError as e:
    print(f"✗ Failed to import SE3Transformer: {e}")

print("\n" + "="*50)
print("Detailed test completed.")
print("="*50)





