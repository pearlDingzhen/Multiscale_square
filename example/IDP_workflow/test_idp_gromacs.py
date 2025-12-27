#!/usr/bin/env python3
"""
Test script for IDP GROMACS runner in example/IDP_workflow
"""

import os
import sys
import yaml
from pathlib import Path

from multiscale2.gromacs_equilibrate import IDPGromacsEquilibrator

def test_idp_gromacs_runner():
    """Test the IDP GROMACS runner with FUS_LC system."""
    
    print("="*80)
    print("Testing IDP GROMACS Runner - FUS_LC System")
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
    gromacs_working_dir = "output_gromacs"
    
    # Check if OpenMM output exists
    if not os.path.exists(openmm_output_dir):
        print(f"✗ OpenMM output directory not found: {openmm_output_dir}")
        print("Please run the OpenMM optimization first:")
        print("  python run_openmm.py")
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
        # Initialize IDP GROMACS equilibrator
        print("Initializing IDP GROMACS equilibrator...")
        runner = IDPGromacsEquilibrator(config, gromacs_working_dir)
        
        # Test Step 1: Get OpenMM outputs
        print("\n" + "="*60)
        print("Testing Step 1: Get OpenMM Outputs")
        print("="*60)
        conf_gro_path, pace_top_path = runner.get_openmm_outputs(openmm_output_dir)
        print(f"✓ Step 1 completed successfully!")
        
        # Test Step 2: Create MDP files
        print("\n" + "="*60)
        print("Testing Step 2: Create MDP Files")
        print("="*60)
        runner._create_em_mdps()
        print(f"✓ Step 2 completed successfully!")
        
        # Check created files
        em_dir = runner.em_dir
        print(f"\nFiles in energy minimization directory:")
        for file in em_dir.iterdir():
            if file.is_file():
                print(f"  - {file.name} ({file.stat().st_size} bytes)")
        
        # Test Step 3: Create equilibration MDP files
        print("\n" + "="*60)
        print("Testing Step 3: Create Equilibration MDP Files")
        print("="*60)
        runner._create_equilibration_mdps()
        print(f"✓ Step 3 completed successfully!")
        
        # Check created files
        eq_dir = runner.eq_dir
        print(f"\nFiles in equilibration directory:")
        for file in eq_dir.iterdir():
            if file.is_file():
                print(f"  - {file.name} ({file.stat().st_size} bytes)")
        
        print(f"\n✓ All tests passed!")
        print(f"✓ IDP GROMACS equilibrator is ready for energy minimization and equilibration!")
        print(f"\nConfiguration loaded:")
        print(f"  - Temperature: {runner.temperature} K")
        print(f"  - Pressure: {runner.pressure} bar")
        print(f"  - CPUs: {runner.num_cpus}")
        print(f"  - NVT steps: {runner.nvt_steps}")
        print(f"  - NPT steps: {runner.npt_steps}")
        print(f"\nNext step: Run complete equilibration workflow")
        print(f"Command would be:")
        print(f"  runner.run_idp_equilibration_workflow('output_openmm/openmm')")
        
        return True
        
    except Exception as e:
        print(f"✗ Error testing IDP GROMACS runner: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main test function."""
    success = test_idp_gromacs_runner()
    
    if success:
        print("\n" + "="*80)
        print("✓ All tests passed!")
        print("✓ IDP GROMACS runner is working correctly!")
        print("="*80)
    else:
        print("\n" + "="*80)
        print("✗ Tests failed!")
        print("="*80)
        sys.exit(1)

if __name__ == "__main__":
    main()
