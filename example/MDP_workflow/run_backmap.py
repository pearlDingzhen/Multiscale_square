#!/usr/bin/env python3
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

# Import required modules for backmapping functions
import multiscale2
sys.path.append(os.path.dirname(multiscale2.__path__[0]))
from multiscale2.backmap import (
    Backmapper, get_latest_pdb, get_protein_dimensions,
    check_protein_fits_in_box, write_pdb_with_bfactors,
    center_condensate, add_chain_ids, remove_ot2_atoms, add_ter_records,
    write_cryst1_only, read_cryst1_dims
)

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

        # Write CRYST1 with new dimensions (convert from nm to Angstroms)
        final_dims_A = new_dims * 10.0  # Convert nm to Angstroms
        cryst1_fixed_path = os.path.join(output_dir, "cryst1_fixed.pdb")
        write_cryst1_only(processed_pdb_path, cryst1_fixed_path, final_dims_A)
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

    # Remove OT2 atoms from backmapped structure
    cleaned_aa_pdb = os.path.join(output_dir, "backmapped_aa_cleaned.pdb")
    remove_ot2_atoms(output_aa_pdb, cleaned_aa_pdb)
    output_aa_pdb = cleaned_aa_pdb  # Use cleaned version for subsequent processing

    # Ensure AA PDB has correct CRYST1 and chain IDs
    print("Finalizing AA structure...")
    if BOX_RESIZE.get('enabled', False):
        # Convert from nm to Angstroms for final dimensions
        final_dims_used = new_dims * 10.0
    else:
        u_src = mda.Universe(source_pdb_path)
        a, b, c = u_src.dimensions[:3]
        final_dims_used = np.array([a, b, c], dtype=float)

    aa_box_pdb = os.path.join(output_dir, "backmapped_aa_box.pdb")
    write_cryst1_only(output_aa_pdb, aa_box_pdb, final_dims_used)

    if os.path.exists(components_path):
        aa_final_pdb = os.path.join(output_dir, "backmapped_aa_final.pdb")
        add_chain_ids(aa_box_pdb, aa_final_pdb, components_path)
        
        # Add TER records to final structure
        aa_final_with_ter = os.path.join(output_dir, "backmapped_aa_final_ter.pdb")
        add_ter_records(aa_final_pdb, aa_final_with_ter)
        final_output = aa_final_with_ter
    else:
        # Add TER records even if no chain IDs
        aa_final_with_ter = os.path.join(output_dir, "backmapped_aa_final_ter.pdb")
        add_ter_records(aa_box_pdb, aa_final_with_ter)
        final_output = aa_final_with_ter

    print("="*60)
    print("Backmapping completed successfully!")
    print(f"Final structure: {final_output}")
    print("="*60)

if __name__ == "__main__":
    run_backmapping()
