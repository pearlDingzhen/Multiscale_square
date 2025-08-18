#!/usr/bin/env python3
"""
Workflow Generator for Multiscale Simulation
Automatically generates run scripts for each stage based on YAML configuration.
"""

import os
import yaml
import glob
from pathlib import Path

class WorkflowGenerator:
    """Generates run scripts for multiscale simulation workflow."""
    
    def __init__(self, config_path):
        """Initialize with configuration file path."""
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.config_path = config_path
        self.config_dir = os.path.dirname(os.path.abspath(config_path))
        
        # Global task type (IDP or MDP). Prefer top-level, fallback to cg_calvados.task_type, default IDP
        self.task_type = (
            self.config.get('task_type')
            or self.config.get('cg_calvados', {}).get('task_type')
            or 'IDP'
        )
        
    def generate_all_scripts(self, output_dir="."):
        """Generate all run scripts for the workflow."""
        print("Generating workflow scripts...")
        
        # Generate Stage 1: CALVADOS simulation
        self.generate_calvados_script(output_dir)
        
        # Generate Stage 2: Backmapping
        self.generate_backmap_script(output_dir)
        
        print("Workflow scripts generated successfully!")
        
    def generate_calvados_script(self, output_dir):
        """Generate run_calvados.py script."""
        script_content = '''#!/usr/bin/env python3
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
    response = input("\\nRun CALVADOS simulation? (y/n): ")
    if response.lower() == 'y':
        print("\\nStarting CALVADOS simulation...")
        
        # Step 1: Run prepare.py
        print("\\nStep 1: Running prepare.py...")
        result1 = subprocess.run(commands[0], shell=True)
        
        if result1.returncode == 0:
            print("✓ prepare.py completed successfully!")
            
            # Step 2: Run the actual simulation
            print("\\nStep 2: Running CALVADOS simulation...")
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
        print("\\nManual commands:")
        for i, cmd in enumerate(commands, 1):
            print(f"  Step {i}: {cmd}")

if __name__ == "__main__":
    run_calvados_simulation()
'''
        
        script_path = os.path.join(output_dir, "run_calvados.py")
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        print(f"Generated: {script_path}")
        
    def generate_backmap_script(self, output_dir):
        """Generate run_backmap.py script."""
        script_content = '''#!/usr/bin/env python3
"""
Stage 2: Backmapping CG to AA
Auto-generated script for backmapping coarse-grained structure to all-atom.
"""

import os
import sys
import yaml
import glob
import numpy as np
import MDAnalysis as mda

# ============================================================================
# USER CONFIGURABLE VARIABLES
# ============================================================================
# These variables can be modified by the user if needed

# Input directory from Stage 1 (leave None to auto-detect)
INPUT_DIR = None

# CG PDB file to backmap (leave "LATEST" to auto-detect, or specify filename)
CG_PDB_FILE = "LATEST"

# Output directory for backmapping (default)
OUTPUT_DIR = "output_backmap"

# Box resize configuration (set enabled=False to disable)
BOX_RESIZE = {
    "enabled": True,
    "new_box_dimensions": [23.0, 17.0, 47.0]  # nm
}

# Note: Tool configuration is automatically determined based on task type:
# - IDP: CA bead + fix_sidechains=True
# - MDP: RES bead + fix_sidechains=False

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

def auto_detect_input_dir():
    """Auto-detect input directory from Stage 1 output."""
    config = load_config()
    protein_config = config.get('protein', {})
    protein_name = protein_config.get('name', 'protein')
    
    # Default CALVADOS output directory: output_calvados/{protein_name}_1
    calvados_dir = os.path.join("output_calvados", f"{protein_name}_1")
    if os.path.exists(calvados_dir):
        return calvados_dir
    
    print(f"Default CALVADOS directory not found: {calvados_dir}")
    return None

def find_suitable_frame_from_trajectory(topology_pdb, trajectory_dcd, box_dims, output_dir, search_last_percent=10.0):
    """Search trajectory for a frame that fits in the specified box."""
    print("Searching trajectory for suitable frame...")
    
    try:
        u = mda.Universe(topology_pdb, trajectory_dcd)
    except Exception as e:
        print(f"Error loading trajectory: {e}")
        return None

    n_frames = len(u.trajectory)
    start_frame = int(n_frames * (1 - search_last_percent / 100.0))
    
    for i in range(n_frames - 1, start_frame - 1, -1):
        u.trajectory[i]
        protein_dims = get_protein_dimensions(u)      # Angstroms
        protein_dims_nm = protein_dims / 10.0         # Convert to nm
        if check_protein_fits_in_box(protein_dims_nm, box_dims):
            suitable_frame_pdb = os.path.join(output_dir, f"suitable_frame_{i}.pdb")
            
            # Set box dimensions directly in nm
            u.dimensions = list(box_dims) + [90, 90, 90]
            
            # Write PDB with correct dimensions
            write_pdb_with_bfactors(u, suitable_frame_pdb)

            print(f"Found suitable frame {i}")
            return suitable_frame_pdb
            
    print(f"No suitable frame found in last {search_last_percent}% of trajectory")
    return None

def run_backmapping():
    """Execute backmapping process."""
    print("="*60)
    print("Stage 2: Backmapping CG to AA")
    print("="*60)
    
    # Load configuration
    config = load_config()
    backmap_config = config.get('backmapping', {})
    
    # Use default output directory
    output_dir = "output_backmap"
    
    # Determine input directory
    input_dir = INPUT_DIR or auto_detect_input_dir()
    if not input_dir:
        print("Error: Could not determine input directory")
        print("Please set INPUT_DIR variable or ensure Stage 1 completed successfully")
        sys.exit(1)
    
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Import required modules
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from multiscale2.backmap import (
        Backmapper, get_latest_pdb, get_protein_dimensions,
        check_protein_fits_in_box, write_pdb_with_bfactors,
        center_condensate, add_chain_ids, enforce_xy_keep_z, 
        write_cryst1_only, read_cryst1_dims
    )
    
    # Locate input CG PDB
    if CG_PDB_FILE == "LATEST":
        source_pdb_path = get_latest_pdb(input_dir)
        if not source_pdb_path:
            print(f"Error: No PDB files found in '{input_dir}'")
            sys.exit(1)
        print(f"Using latest PDB: {os.path.basename(source_pdb_path)}")
    else:
        source_pdb_path = os.path.join(input_dir, CG_PDB_FILE)
        if not os.path.exists(source_pdb_path):
            print(f"Error: PDB file not found: {source_pdb_path}")
            sys.exit(1)
    
    # Pre-processing: center and add chain IDs
    centered_pdb = center_condensate(
        source_pdb_path,
        os.path.join(output_dir, "centered_cg.pdb")
    )

    components_path = os.path.join(input_dir, 'components.yaml')
    if os.path.exists(components_path):
        chained_pdb = add_chain_ids(
            centered_pdb,
            os.path.join(output_dir, "chained_cg.pdb"),
            components_path
        )
        processed_pdb_path = chained_pdb
    else:
        processed_pdb_path = centered_pdb
        print("Warning: components.yaml not found, skipping chain ID assignment")
    
    # Box resize if enabled
    if BOX_RESIZE.get('enabled', False):
        print("Validating structure against new box dimensions...")
        # Prefer values from config.yaml if present
        cfg_dims = backmap_config.get('box_resize', {}).get('new_box_dimensions') if backmap_config else None
        new_dims = np.array(cfg_dims if cfg_dims is not None else BOX_RESIZE.get('new_box_dimensions'))
        
        u = mda.Universe(processed_pdb_path)
        protein_dims = get_protein_dimensions(u)  # Angstroms
        protein_dims_nm = protein_dims / 10.0      # Convert to nm
        print(f"Protein dimensions (Angstroms): {protein_dims}")
        print(f"Protein dimensions (nm): {protein_dims_nm}")
        print(f"Target box dimensions (nm): {new_dims}")

        if not check_protein_fits_in_box(protein_dims_nm, new_dims):
            print("Protein dimensions exceed box size. Searching trajectory...")
            
            topology_file = os.path.join(input_dir, 'top.pdb')
            dcd_files = glob.glob(os.path.join(input_dir, '*.dcd'))
            
            if os.path.exists(topology_file) and dcd_files:
                trajectory_file = max(dcd_files, key=os.path.getctime)
                processed_pdb_path = find_suitable_frame_from_trajectory(
                    topology_file, trajectory_file, new_dims, output_dir
                )
                if not processed_pdb_path:
                    print("No suitable frame found in trajectory.")
                    print("Please modify BOX_RESIZE['new_box_dimensions'] in the script header, or set CG_PDB_FILE to a specific PDB filename.")
                    sys.exit(1)
            else:
                print("Error: Topology or trajectory files not found")
                print("Please modify BOX_RESIZE['new_box_dimensions'] in the script header, or set CG_PDB_FILE to a specific PDB filename.")
                sys.exit(1)
        
        # Resize box
        resized_pdb_path = os.path.join(output_dir, "resized_cg.pdb")
        final_universe = mda.Universe(processed_pdb_path)
        # Convert from nm to nm (no conversion needed since new_dims is already in nm)
        final_universe.dimensions = list(new_dims) + [90, 90, 90]
        write_pdb_with_bfactors(final_universe, resized_pdb_path)
        processed_pdb_path = resized_pdb_path

        # Write CRYST1 with correct dimensions
        final_dims = enforce_xy_keep_z(source_pdb_path, new_dims)
        cryst1_fixed_path = os.path.join(output_dir, "cryst1_fixed.pdb")
        write_cryst1_only(processed_pdb_path, cryst1_fixed_path, final_dims)
        processed_pdb_path = cryst1_fixed_path
    
    # Determine tool configuration based on task type
    config = load_config()
    # Prefer global task_type, fallback to legacy location
    task_type = config.get('task_type', config.get('backmapping', {}).get('task_type', 'IDP'))
    
    if task_type == 'IDP':
        tool_config = {
            "tool": "cg2all",
            "cg_bead_name": "CA",
            "fix_sidechains": True
        }
    elif task_type == 'MDP':
        tool_config = {
            "tool": "cg2all",
            "cg_bead_name": "RES",
            "fix_sidechains": False
        }
    else:
        print(f"Warning: Unknown task type '{task_type}', using IDP defaults")
        tool_config = {
            "tool": "cg2all",
            "cg_bead_name": "CA",
            "fix_sidechains": True
        }
    
    print(f"Using tool configuration for {task_type}: {tool_config}")
    
    # Perform backmapping
    output_aa_pdb = os.path.join(output_dir, "backmapped_aa.pdb")
    
    try:
        mapper = Backmapper(tool_config)
        mapper.reconstruct(
            input_cg_pdb=processed_pdb_path,
            output_aa_pdb=output_aa_pdb
        )
    except Exception as e:
        print(f"Error during backmapping: {e}")
        sys.exit(1)

    # Ensure AA PDB has correct CRYST1 and chain IDs
    print("Finalizing AA structure...")
    if BOX_RESIZE.get('enabled', False):
        final_dims_used = enforce_xy_keep_z(source_pdb_path, new_dims)
    else:
        u_src = mda.Universe(source_pdb_path)
        a, b, c = u_src.dimensions[:3]
        final_dims_used = np.array([a, b, c], dtype=float)

    aa_box_pdb = os.path.join(output_dir, "backmapped_aa_box.pdb")
    write_cryst1_only(output_aa_pdb, aa_box_pdb, final_dims_used)

    if os.path.exists(components_path):
        aa_final_pdb = os.path.join(output_dir, "backmapped_aa_final.pdb")
        add_chain_ids(aa_box_pdb, aa_final_pdb, components_path)
        final_output = aa_final_pdb
    else:
        final_output = aa_box_pdb

    print("="*60)
    print("Backmapping completed successfully!")
    print(f"Final structure: {final_output}")
    print("="*60)

if __name__ == "__main__":
    run_backmapping()
'''
        
        script_path = os.path.join(output_dir, "run_backmap.py")
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        print(f"Generated: {script_path}")

def main():
    """Main function to generate workflow scripts."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate multiscale workflow scripts")
    parser.add_argument("config", help="Path to configuration YAML file")
    parser.add_argument("--output", "-o", default=".", help="Output directory for scripts")
    
    args = parser.parse_args()
    
    generator = WorkflowGenerator(args.config)
    generator.generate_all_scripts(args.output)

if __name__ == "__main__":
    main()
