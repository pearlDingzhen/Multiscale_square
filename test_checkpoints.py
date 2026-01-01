#!/usr/bin/env python
"""
Test script to verify checkpoint loading from the provided model directory.
"""

import sys
sys.path.insert(0, '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square')

from multiscale2.extern.ms2_cg2all.lib.libconfig import MODEL_HOME
import os

print("="*60)
print("Testing Checkpoint Loading")
print("="*60)

# List available checkpoints
print(f"\nMODEL_HOME: {MODEL_HOME}")
print(f"Directory exists: {MODEL_HOME.exists()}")
print(f"\nAvailable checkpoints:")

available_models = {}
for f in sorted(MODEL_HOME.glob("*.ckpt")):
    file_size = f.stat().st_size
    model_name = f.name.replace("-FIX.ckpt", "").replace(".ckpt", "")
    
    if model_name not in available_models:
        available_models[model_name] = {}
    
    if "-FIX" in f.name:
        available_models[model_name]["FIX"] = file_size
    else:
        available_models[model_name]["normal"] = file_size
    
    size_str = f"{file_size / (1024*1024):.1f} MB" if file_size > 1024*1024 else f"{file_size / 1024:.1f} KB"
    print(f"  {f.name}: {size_str}")

print(f"\n{len(available_models)} unique models found")

# Test loading a checkpoint
print("\n" + "="*60)
print("Testing checkpoint loading...")
print("="*60)

import torch

test_models = ["CalphaBasedModel", "ResidueBasedModel", "Martini", "Martini3"]
for model_name in test_models:
    for variant in ["normal", "FIX"]:
        if variant == "FIX":
            ckpt_name = f"{model_name}-FIX.ckpt"
        else:
            ckpt_name = f"{model_name}.ckpt"
        
        ckpt_path = MODEL_HOME / ckpt_name
        if not ckpt_path.exists():
            print(f"  {ckpt_name}: NOT FOUND (skipping)")
            continue
        
        try:
            file_size = ckpt_path.stat().st_size
            if file_size < 1000:
                print(f"  {ckpt_name}: FILE TOO SMALL ({file_size} bytes) - skipping")
                continue
                
            print(f"  {ckpt_name}: Loading...", end=" ")
            ckpt = torch.load(ckpt_path, map_location="cpu")
            
            # Verify checkpoint contents
            if "hyper_parameters" in ckpt:
                cg_model = ckpt["hyper_parameters"].get("cg_model", "unknown")
                print(f"OK (cg_model={cg_model})")
            else:
                print("WARNING: No hyper_parameters found")
        except Exception as e:
            print(f"ERROR: {e}")

print("\n" + "="*60)
print("Checkpoint test completed!")
print("="*60)





