#!/usr/bin/env python3
"""
Example script demonstrating how to use the workflow generator
to create run scripts for multiscale simulation.
"""

import os
import sys

# Add the multiscale2 module to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from multiscale2.workflow_generator import WorkflowGenerator

def main():
    """Generate workflow scripts for the IDP example."""
    print("="*60)
    print("Multiscale Workflow Generator Example")
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
    print("\nNext steps:")
    print("1. Review and modify the generated scripts if needed")
    print("2. Run Stage 1: python run_calvados.py")
    print("3. Run Stage 2: python run_backmap.py")
    print("\nNote: Each script contains user-configurable variables at the top")
    print("that can be modified before running.")

if __name__ == "__main__":
    main()
