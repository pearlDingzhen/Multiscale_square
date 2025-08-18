#!/usr/bin/env python3
"""
Stage 1: CALVADOS Coarse-Grained Simulation
Auto-generated script for running CALVADOS simulation.
"""

import os
import sys
import subprocess
import yaml

# ============================================================================
# USER CONFIGURABLE VARIABLES
# ============================================================================
# These variables can be modified by the user if needed

# Output directory for CALVADOS simulation
OUTPUT_DIR = "output_calvados"

# GPU device ID (0 for first GPU)
GPU_ID = 0

# Replica number
REPLICA = 1

# Override protein name from config (leave None to use config value)
PROTEIN_NAME_OVERRIDE = None

# Override number of protein molecules (leave None to use config value)
PROTEIN_NMOL_OVERRIDE = None

# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

def load_config():
    """Load configuration from YAML file."""
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def run_calvados_simulation():
    """Execute CALVADOS simulation."""
    print("="*60)
    print("Stage 1: CALVADOS Coarse-Grained Simulation")
    print("="*60)
    
    # Load configuration
    config = load_config()
    calvados_config = config.get('cg_calvados', {})
    protein_config = config.get('protein', {})
    
    # Get protein settings (with override support)
    protein_name = PROTEIN_NAME_OVERRIDE or protein_config.get('name', 'protein')
    protein_nmol = PROTEIN_NMOL_OVERRIDE or protein_config.get('nmol', 20)
    
    print(f"Protein name: {protein_name}")
    print(f"Number of molecules: {protein_nmol}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"GPU ID: {GPU_ID}")
    print(f"Replica: {REPLICA}")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Import and initialize generator
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from multiscale2.calvados_generator import CalvadosGenerator
    
    generator = CalvadosGenerator("config.yaml")
    
    # Override protein settings if specified
    if PROTEIN_NAME_OVERRIDE:
        generator.protein_name = PROTEIN_NAME_OVERRIDE
    if PROTEIN_NMOL_OVERRIDE:
        generator.protein_nmol = PROTEIN_NMOL_OVERRIDE
    
    # Generate CALVADOS commands
    commands = generator.generate_and_run(OUTPUT_DIR, protein_name, gpu_id=GPU_ID, replica=REPLICA)
    
    print(f"Generated commands:")
    for i, cmd in enumerate(commands, 1):
        print(f"  Step {i}: {cmd}")
    
    # Ask user if they want to run the simulation
    response = input("\nRun CALVADOS simulation? (y/n): ")
    if response.lower() == 'y':
        print("\nStarting CALVADOS simulation...")
        
        # Step 1: Run prepare.py
        print("\nStep 1: Running prepare.py...")
        result1 = subprocess.run(commands[0], shell=True)
        
        if result1.returncode == 0:
            print("✓ prepare.py completed successfully!")
            
            # Step 2: Run the actual simulation
            print("\nStep 2: Running CALVADOS simulation...")
            result2 = subprocess.run(commands[1], shell=True)
            
            if result2.returncode == 0:
                print("✓ CALVADOS simulation completed successfully!")
                
                # Check for output files
                calvados_dir = os.path.join(OUTPUT_DIR, f"{protein_name}_{REPLICA}")
                if os.path.exists(calvados_dir):
                    print(f"Output directory: {calvados_dir}")
                    files = os.listdir(calvados_dir)
                    print(f"Generated {len(files)} files")
                    
                    # Look for trajectory and structure files
                    trajectory_files = [f for f in files if f.endswith('.dcd')]
                    structure_files = [f for f in files if f.endswith('.pdb')]
                    
                    if trajectory_files:
                        print(f"Trajectory files: {trajectory_files}")
                    if structure_files:
                        print(f"Structure files: {structure_files}")
            else:
                print("❌ CALVADOS simulation failed!")
                sys.exit(1)
        else:
            print("❌ prepare.py failed!")
            sys.exit(1)
    else:
        print("Simulation skipped. You can run it manually later.")
        print("\nManual commands:")
        for i, cmd in enumerate(commands, 1):
            print(f"  Step {i}: {cmd}")

if __name__ == "__main__":
    run_calvados_simulation()
