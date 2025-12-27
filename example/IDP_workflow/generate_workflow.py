#!/usr/bin/env python3
"""
Example script demonstrating how to use the workflow generator
to create run scripts for multiscale simulation.
"""

import os
import sys

from multiscale2.workflow_generator import WorkflowGenerator

def main():
    """Generate workflow scripts for the IDP example."""
    print("="*60)
    print("Multiscale Workflow Generator Example (IDP)")
    print("="*60)
    
    # Configuration file path
    config_path = "config.yaml"
    
    if not os.path.exists(config_path):
        print(f"Error: Configuration file not found: {config_path}")
        print("Please ensure you are in the example/IDP directory with config.yaml")
        return
    
    print(f"Using configuration file: {config_path}")
    
    # Create workflow generator
    generator = WorkflowGenerator(config_path)
    
    # Generate scripts in current directory
    output_dir = "."
    generator.generate_all_scripts(output_dir)
    
    print("\n" + "="*60)
    print("Workflow Generation Complete!")
    print("="*60)
    print("\nGenerated files:")
    print("  - run_calvados.py    (Stage 1: CALVADOS simulation)")
    print("  - run_backmap.py     (Stage 2: Backmapping)")
    print("  - run_openmm.py      (Stage 3: OpenMM optimization)")
    print("\nNext steps:")
    print("1. Review and modify the generated scripts if needed")
    print("2. Run Stage 1: python run_calvados.py")
    print("3. Run Stage 2: python run_backmap.py")
    print("4. Run Stage 3: python run_openmm.py")
    print("\nNote: Each script contains user-configurable variables at the top")
    print("that can be modified before running.")
    
    # Test the generated run_openmm.py script
    print("\n" + "="*60)
    print("Testing run_openmm.py script...")
    print("="*60)
    
    test_openmm_script()

def test_openmm_script():
    """Test the generated run_openmm.py script for basic functionality."""
    import subprocess
    import yaml
    
    # Load config to get task type
    with open("config.yaml", 'r') as f:
        config = yaml.safe_load(f)
    
    task_type = config.get('task_type', 'IDP')
    print(f"Task type: {task_type}")
    
    # Check if run_openmm.py exists
    if not os.path.exists("run_openmm.py"):
        print("❌ run_openmm.py not found!")
        return
    
    print("✓ run_openmm.py found")
    
    # Test basic syntax
    try:
        result = subprocess.run([sys.executable, "-m", "py_compile", "run_openmm.py"], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ run_openmm.py syntax is valid")
        else:
            print(f"❌ run_openmm.py has syntax errors:")
            print(result.stderr)
            return
    except Exception as e:
        print(f"❌ Error testing run_openmm.py: {e}")
        return
    
    # Test imports
    try:
        # Test if the script can import required modules
        test_code = '''
import sys
import os
try:
    from multiscale2 import PACE2openmm
    print("✓ PACE2openmm import successful")
except ImportError as e:
    print(f"❌ PACE2openmm import failed: {e}")
    sys.exit(1)
'''
        result = subprocess.run([sys.executable, "-c", test_code], 
                              capture_output=True, text=True)
        print(result.stdout.strip())
        if result.returncode != 0:
            print(result.stderr.strip())
            return
    except Exception as e:
        print(f"❌ Error testing imports: {e}")
        return
    
    # Test configuration loading
    try:
        test_code = '''
import yaml
with open("config.yaml", 'r') as f:
    config = yaml.safe_load(f)
task_type = config.get('task_type', 'IDP')
print(f"✓ Config loaded successfully, task_type: {task_type}")
'''
        result = subprocess.run([sys.executable, "-c", test_code], 
                              capture_output=True, text=True)
        print(result.stdout.strip())
        if result.returncode != 0:
            print(result.stderr.strip())
            return
    except Exception as e:
        print(f"❌ Error testing config loading: {e}")
        return
    
    print("\n" + "="*60)
    print("✅ run_openmm.py test completed successfully!")
    print("="*60)
    print("\nThe script is ready to use. You can now:")
    print("1. Complete Stage 1 and Stage 2 first")
    print("2. Run: python run_openmm.py")
    print("\nThe script will automatically:")
    print(f"  - Generate PACE topology for {task_type} task")
    print(f"  - Process structure with pdb2gmx")
    print(f"  - Run OpenMM optimization {'with CA restraints' if task_type == 'IDP' else 'without CA restraints'}")

if __name__ == "__main__":
    main()
