#!/usr/bin/env python3
"""
Test 1: COCOMO Simulation with config_idp.yaml

This test runs a complete COCOMO simulation using the IDP configuration.
The system consists of 100 TDP43_CTD chains in a cubic box.

Note: This test requires GPU for reasonable performance.
Set gpu_id=-1 to use CPU (slower but works without GPU).
"""

import os
import sys
import time
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

from multiscale2.src import CGSimulationConfig, CGSimulator


def main():
    """Run COCOMO simulation with config_idp.yaml."""
    print("\n" + "=" * 60)
    print("COCOMO2 Creator - Test 1: Full Simulation")
    print("=" * 60 + "\n")
    
    # Configuration
    test_dir = Path(__file__).parent
    config_file = test_dir / "config_idp.yaml"
    output_dir = test_dir / "output_test1"
    
    # Simulation parameters
    gpu_id = 0  # Set to -1 for CPU, 0+ for GPU
    preequil_steps = 10000  # Pre-equilibration steps
    simulation_steps = 500000  # Main simulation steps (increased for better sampling)
    
    print(f"Configuration file: {config_file}")
    print(f"Output directory: {output_dir}")
    print(f"GPU ID: {gpu_id}")
    print(f"Pre-equilibration steps: {preequil_steps}")
    print(f"Simulation steps: {simulation_steps}")
    print()
    
    # Check if config file exists
    if not config_file.exists():
        print(f"[ERROR] Configuration file not found: {config_file}")
        return 1
    
    # Load configuration
    print("Loading configuration...")
    try:
        config = CGSimulationConfig.from_yaml(str(config_file))
        print(f"  System name: {config.system_name}")
        print(f"  Box size: {config.box} nm")
        print(f"  Temperature: {config.temperature} K")
        print(f"  Components: {len(config.components)}")
        for comp in config.components:
            print(f"    - {comp.name}: {comp.type.value}, nmol={comp.nmol}")
        print()
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1
    
    # Create simulator
    print("Creating simulator...")
    try:
        sim = CGSimulator(config)
        print("  Simulator created successfully")
        print()
    except Exception as e:
        print(f"[ERROR] Failed to create simulator: {e}")
        return 1
    
    # Setup output directory
    print("Setting up output directory...")
    try:
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        sim.setup(str(output_dir), overwrite=True)
        print(f"  Output directory: {output_dir}")
        print()
    except Exception as e:
        print(f"[ERROR] Failed to setup simulator: {e}")
        return 1
    
    # Print topology info
    print("Topology Information:")
    print(f"  Global sequence length: {len(sim.get_global_sequence())}")
    print(f"  Number of chains: {max(sim.get_chain_ids())}")
    print(f"  Folded domains: {sum(sim.get_folded_domains())}")
    print()
    
    # Run COCOMO simulation
    print("=" * 60)
    print("Starting COCOMO Simulation")
    print("=" * 60)
    
    start_time = time.time()

    try:
        # Set simulation steps before running
        original_steps = sim.config.simulation.steps
        sim.config.simulation.steps = simulation_steps
        print(f"  Simulation steps: {simulation_steps} (config: {original_steps})")

        result = sim.run_cocomo(
            gpu_id=gpu_id,
            preequil_steps=preequil_steps,
        )

        # Verify output files
        if result.success:
            print()
            print("Verifying output files...")

            output_files = [
                'trajectory.xtc',
                'final.pdb',
                'system.xml',
                'simulation.log',
            ]

            all_exist = True
            for f in output_files:
                filepath = output_dir / f
                if filepath.exists():
                    size = filepath.stat().st_size
                    print(f"  [EXISTS] {f} ({size:,} bytes)")
                else:
                    print(f"  [MISSING] {f}")
                    all_exist = False

            if all_exist:
                print("\n[SUCCESS] All output files generated!")
            else:
                print("\n[WARNING] Some output files are missing")
        
        elapsed = time.time() - start_time
        print(f"\nTotal simulation time: {elapsed:.1f} seconds")
        
        if result.success:
            print("\n[SUCCESS] COCOMO simulation completed successfully!")
            return 0
        else:
            print(f"\n[FAILURE] Simulation failed: {result.errors}")
            return 1
            
    except KeyboardInterrupt:
        print("\n[WARNING] Simulation interrupted by user")
        return 130
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n[ERROR] Simulation failed after {elapsed:.1f} seconds")
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    exit(main())

