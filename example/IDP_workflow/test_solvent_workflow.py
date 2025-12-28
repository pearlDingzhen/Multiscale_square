#!/usr/bin/env python3
"""
Test script for the solvent preparation workflow
"""

import os
import sys
import yaml
from pathlib import Path

from run_solvent import SolventPreparator

def test_solvent_workflow():
    """Test the solvent preparation workflow with IDP system."""
    
    print("="*80)
    print("Testing Solvent Preparation Workflow - FUS_LC System")
    print("="*80)
    
    # Load configuration
    config_path = "config.yaml"
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    print(f"Testing with protein: {config['protein']['name']}")
    print(f"Number of molecules: {config['protein']['nmol']}")
    print(f"Task type: {config['task_type']}")
    print()
    
    # Set up paths
    openmm_output_dir = "output_openmm/openmm"
    working_dir = "test_solvent_output"
    
    # Check if OpenMM output exists
    if not os.path.exists(openmm_output_dir):
        print(f"✗ OpenMM output directory not found: {openmm_output_dir}")
        print("Please run the OpenMM optimization first:")
        print("  cd example/IDP_workflow && python run_openmm.py")
        return False
    
    # Check for required files
    conf_gro = os.path.join(openmm_output_dir, "conf.gro")
    pace_top = os.path.join(openmm_output_dir, "PACE.top")
    
    if not os.path.exists(conf_gro):
        print(f"✗ conf.gro not found: {conf_gro}")
        return False
    
    if not os.path.exists(pace_top):
        print(f"✗ PACE.top not found: {pace_top}")
        return False
    
    print(f"✓ Found OpenMM outputs:")
    print(f"  - {conf_gro}")
    print(f"  - {pace_top}")
    print()
    
    try:
        # Initialize solvent preparator
        print("Initializing solvent preparator...")
        preparator = SolventPreparator(config, working_dir)
        
        # Test Step 1: Vacuum optimization setup
        print("\n" + "="*60)
        print("Testing Step 1: Vacuum Optimization Setup")
        print("="*60)
        preparator._create_vacuum_em_mdps()
        print(f"✓ Step 1 setup completed successfully!")
        
        # Test Step 2: Solvent system setup
        print("\n" + "="*60)
        print("Testing Step 2: Solvent System Setup")
        print("="*60)
        preparator._create_solvent_em_mdps()
        print(f"✓ Step 2 setup completed successfully!")
        
        # Check created files
        vacuum_dir = preparator.vacuum_dir
        solvent_dir = preparator.solvent_dir
        
        print(f"\nFiles in vacuum optimization directory:")
        for file in vacuum_dir.iterdir():
            if file.is_file():
                print(f"  - {file.name} ({file.stat().st_size} bytes)")
        
        print(f"\nFiles in explicit solvent directory:")
        for file in solvent_dir.iterdir():
            if file.is_file():
                print(f"  - {file.name} ({file.stat().st_size} bytes)")
        
        print(f"\n✓ All tests passed!")
        print(f"✓ Solvent preparator is ready for full workflow!")
        print(f"\nConfiguration loaded:")
        print(f"  - Temperature: {preparator.temperature} K")
        print(f"  - Pressure: {preparator.pressure} bar")
        print(f"  - CPUs: {preparator.num_cpus}")
        print(f"  - Ions concentration: {preparator.ions_concentration} M")
        print(f"  - Anti-freeze water: {preparator.use_anti_freeze_water}")
        print(f"  - Anti-freeze ratio: {preparator.anti_free_water_ratio}")
        print(f"\nNext step: Run complete solvent preparation workflow")
        print(f"Command would be:")
        print(f"  python run_solvent.py")
        print(f"  # or")
        print(f"  python run_solvent.py {config_path} {openmm_output_dir}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error testing solvent preparator: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main test function."""
    success = test_solvent_workflow()
    
    if success:
        print("\n" + "="*80)
        print("✓ All tests passed!")
        print("✓ Solvent preparation workflow is ready!")
        print("="*80)
    else:
        print("\n" + "="*80)
        print("✗ Tests failed!")
        print("="*80)
        sys.exit(1)

if __name__ == "__main__":
    main()
