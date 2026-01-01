#!/usr/bin/env python
"""
Test script for ms2_cg2all integration.
This script tests the full import chain from multiscale2.extern.ms2_cg2all
"""

import sys
import os

# Add the project root to the path
sys.path.insert(0, '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square')

print("="*60)
print("Testing ms2_cg2all Integration")
print("="*60)

# Test 1: Import the main module
print("\n1. Testing main module import...")
try:
    from multiscale2.extern.ms2_cg2all import convert_cg2all
    print("✓ Successfully imported convert_cg2all from multiscale2.extern.ms2_cg2all")
except ImportError as e:
    print(f"✗ Failed to import convert_cg2all: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Import SE3Transformer
print("\n2. Testing SE3Transformer import...")
try:
    from multiscale2.extern.ms2_se3transformer import Fiber, SE3Transformer
    print("✓ Successfully imported Fiber and SE3Transformer")
except ImportError as e:
    print(f"✗ Failed to import SE3Transformer: {e}")

# Test 3: Test convert_cg2all function signature
print("\n3. Testing convert_cg2all function signature...")
try:
    import inspect
    sig = inspect.signature(convert_cg2all)
    params = list(sig.parameters.keys())
    print(f"✓ Function parameters: {params}")
    required_params = ['in_pdb_fn', 'out_fn']
    for param in required_params:
        if param in params:
            print(f"  ✓ {param} is present")
        else:
            print(f"  ✗ {param} is missing!")
except Exception as e:
    print(f"✗ Failed to check function signature: {e}")

# Test 4: Check that MODEL_HOME is configured correctly
print("\n4. Testing configuration paths...")
try:
    from multiscale2.extern.ms2_cg2all.lib.libconfig import BASE, DATA_HOME, MODEL_HOME
    print(f"✓ BASE: {BASE}")
    print(f"✓ DATA_HOME exists: {DATA_HOME.exists()}")
    print(f"✓ MODEL_HOME: {MODEL_HOME}")
except Exception as e:
    print(f"✗ Failed to check configuration: {e}")

# Test 5: Check data files
print("\n5. Checking data files...")
try:
    from multiscale2.extern.ms2_cg2all.lib.libconfig import DATA_HOME
    data_files = list(DATA_HOME.glob("*.dat")) + list(DATA_HOME.glob("*.pkl")) + list(DATA_HOME.glob("*.top"))
    print(f"✓ Found {len(data_files)} data files in DATA_HOME")
    for f in data_files[:5]:
        print(f"  - {f.name}")
    if len(data_files) > 5:
        print(f"  ... and {len(data_files) - 5} more")
except Exception as e:
    print(f"✗ Failed to check data files: {e}")

print("\n" + "="*60)
print("Integration test completed!")
print("="*60)
