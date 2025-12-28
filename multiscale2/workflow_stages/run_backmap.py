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
from multiscale2.backmap import (
    Backmapper, get_latest_pdb, get_protein_dimensions,
    check_protein_fits_in_box, write_pdb_with_bfactors,
    add_chain_ids, remove_ot2_atoms, add_ter_records,
    write_cryst1_only, read_cryst1_dims
)
import subprocess

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

# Box resize configuration (auto-configured from config.yaml)
BOX_RESIZE = {
    "enabled": False,
    "new_box_dimensions": [8.0, 8.0, 40.0]  # nm
}

# Note: Tool configuration is automatically determined based on task type:
# - IDP: CA bead + fix_sidechains=True
# - MDP: RES bead + fix_sidechains=False

# ============================================================================
# GMX EDITCONF FUNCTIONS
# ============================================================================

def center_condensate_gmx(input_pdb, output_pdb):
    """
    使用gmx editconf居中蛋白质组装体，使用cubic盒子并根据z坐标进行特殊居中
    """
    print(f"--- 使用gmx editconf居中凝聚体: {os.path.basename(input_pdb)} ---")
    
    try:
        # 使用MDAnalysis分析z坐标和盒子尺寸
        u = mda.Universe(input_pdb)
        atoms = u.atoms.copy()
        z_coords = atoms.positions[:, 2]
        sorted_z = np.sort(z_coords)
        center = np.mean(sorted_z)
        
        # z_max是盒子的z方向总长度，不是原子的最大z坐标
        z_max = u.dimensions[2]  # 盒子的z方向长度
        
        z_translate = 0.1 * (0.5 * z_max - center)
        
        print(f"原子Z坐标范围: {sorted_z[0]:.2f} - {sorted_z[-1]:.2f}")
        print(f"盒子Z方向长度: {z_max:.2f}")
        print(f"蛋白质Z中心: {center:.2f}")
        print(f"盒子Z中点: {0.5 * z_max:.2f}")
        print(f"Z平移量: {z_translate:.2f}")
        
        # 使用gmx editconf进行居中操作，使用cubic盒子
        cmd = [
            'gmx', 'editconf',
            '-f', input_pdb,
            '-o', output_pdb,
            '-c',  # 居中分子
            '-translate', '0', '0', str(z_translate),  # 在z方向进行额外平移
            '-bt', 'cubic'  # 设置盒子类型为cubic
        ]
        
        print(f"执行命令: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if os.path.exists(output_pdb):
            print(f"✓ 凝聚体已居中。新PDB保存到: {output_pdb}")
            return output_pdb
        else:
            print(f"❌ 输出文件未生成: {output_pdb}")
            return input_pdb
            
    except subprocess.CalledProcessError as e:
        print(f"❌ gmx editconf居中失败: {e}")
        print(f"错误输出: {e.stderr}")
        return input_pdb
    except FileNotFoundError:
        print("❌ 未找到gmx命令，请确保GROMACS已正确安装并在PATH中")
        return input_pdb
    except Exception as e:
        print(f"❌ 居中过程中发生错误: {e}")
        return input_pdb

def resize_box_gmx(input_pdb, output_pdb, new_box_dims_nm):
    """
    使用gmx editconf调整盒子尺寸，使用cubic盒子，只做resize不做居中
    new_box_dims_nm: 新的盒子尺寸，单位为nm，格式为[x, y, z]
    """
    print(f"--- 使用gmx editconf调整盒子尺寸: {os.path.basename(input_pdb)} ---")
    print(f"新盒子尺寸 (nm): {new_box_dims_nm}")
    
    try:
        # 将nm转换为Angstroms（gmx editconf使用Angstroms）
        new_dims_A = np.array(new_box_dims_nm) * 10.0
        
        # 使用gmx editconf调整盒子尺寸，使用cubic盒子，只做resize不做居中
        cmd = [
            'gmx', 'editconf',
            '-f', input_pdb,
            '-o', output_pdb,
            '-box', str(new_dims_A[0]), str(new_dims_A[1]), str(new_dims_A[2]),
            '-bt', 'cubic'  # 设置盒子类型为cubic
        ]
        
        print(f"执行命令: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if os.path.exists(output_pdb):
            print(f"✓ 盒子尺寸已调整。新PDB保存到: {output_pdb}")
            return output_pdb
        else:
            print(f"❌ 输出文件未生成: {output_pdb}")
            return input_pdb
            
    except subprocess.CalledProcessError as e:
        print(f"❌ gmx editconf调整盒子尺寸失败: {e}")
        print(f"错误输出: {e.stderr}")
        return input_pdb
    except FileNotFoundError:
        print("❌ 未找到gmx命令，请确保GROMACS已正确安装并在PATH中")
        return input_pdb

def check_protein_fits_gmx(pdb_file, box_dims_nm):
    """
    使用MDAnalysis检查蛋白质是否适合新的盒子尺寸
    """
    try:
        u = mda.Universe(pdb_file)
        protein_dims = get_protein_dimensions(u)  # Angstroms
        protein_dims_nm = protein_dims / 10.0     # Convert to nm
        print(f"蛋白质尺寸 (Angstroms): {protein_dims}")
        print(f"蛋白质尺寸 (nm): {protein_dims_nm}")
        print(f"目标盒子尺寸 (nm): {box_dims_nm}")
        
        fits = check_protein_fits_in_box(protein_dims_nm, np.array(box_dims_nm))
        print(f"蛋白质是否适合盒子: {'是' if fits else '否'}")
        return fits
    except Exception as e:
        print(f"❌ 检查蛋白质尺寸失败: {e}")
        return False

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
    
    # Pre-processing: add chain IDs first
    components_path = os.path.join(input_dir, 'components.yaml')
    if os.path.exists(components_path):
        chained_pdb = add_chain_ids(
            source_pdb_path,
            os.path.join(output_dir, "chained_cg.pdb"),
            components_path
        )
        processed_pdb_path = chained_pdb
    else:
        processed_pdb_path = source_pdb_path
        print("Warning: components.yaml not found, skipping chain ID assignment")
    
    # Box resize if enabled
    if BOX_RESIZE.get('enabled', False):
        print("使用gmx editconf验证结构并调整盒子尺寸...")
        # Prefer values from config.yaml if present
        cfg_dims = backmap_config.get('box_resize', {}).get('new_box_dimensions') if backmap_config else None
        new_dims = np.array(cfg_dims if cfg_dims is not None else BOX_RESIZE.get('new_box_dimensions'))
        
        # 使用gmx版本检查蛋白质是否适合盒子
        if not check_protein_fits_gmx(processed_pdb_path, new_dims):
            print("蛋白质尺寸超出盒子大小。搜索轨迹...")
            
            topology_file = os.path.join(input_dir, 'top.pdb')
            dcd_files = glob.glob(os.path.join(input_dir, '*.dcd'))
            
            if os.path.exists(topology_file) and dcd_files:
                trajectory_file = max(dcd_files, key=os.path.getctime)
                processed_pdb_path = find_suitable_frame_from_trajectory(
                    topology_file, trajectory_file, new_dims, output_dir
                )
                if not processed_pdb_path:
                    print("在轨迹中未找到合适的帧。")
                    print("请修改脚本头部的BOX_RESIZE['new_box_dimensions']，或将CG_PDB_FILE设置为特定的PDB文件名。")
                    sys.exit(1)
            else:
                print("错误：未找到拓扑或轨迹文件")
                print("请修改脚本头部的BOX_RESIZE['new_box_dimensions']，或将CG_PDB_FILE设置为特定的PDB文件名。")
                sys.exit(1)
        
        # 第一步：使用gmx editconf调整盒子尺寸
        resized_pdb_path = os.path.join(output_dir, "resized_cg.pdb")
        processed_pdb_path = resize_box_gmx(processed_pdb_path, resized_pdb_path, new_dims)
        
        # 第二步：使用gmx editconf居中
        centered_pdb_path = os.path.join(output_dir, "centered_cg.pdb")
        processed_pdb_path = center_condensate_gmx(processed_pdb_path, centered_pdb_path)
    else:
        # 如果没有resize，只做居中
        centered_pdb_path = os.path.join(output_dir, "centered_cg.pdb")
        processed_pdb_path = center_condensate_gmx(processed_pdb_path, centered_pdb_path)
    
    # Determine tool configuration based on cg_type from config
    config = load_config()
    cg_type = config.get('backmapping', {}).get('cg_type', 'CA')

    # Only CA model needs fix_sidechains, others (martini3, RES, CALVADOS3) don't
    fix_sidechains = (cg_type == 'CA')

    tool_config = {
        "tool": "cg2all",
        "cg_bead_name": cg_type,  # 直接使用 cg_type
        "fix_sidechains": fix_sidechains
    }

    print(f"Using tool configuration: {tool_config}")
    
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
    print("完成AA结构处理...")
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
        
        # Add TER records to structure (save intermediate with OT2 atoms)
        aa_final_with_ter_and_ot2 = os.path.join(output_dir, "backmapped_aa_final_with_OT2.pdb")
        add_ter_records(aa_final_pdb, aa_final_with_ter_and_ot2)
        print(f"Saved intermediate file with OT2 atoms: {aa_final_with_ter_and_ot2}")
        
        # Remove OT2 atoms to create final output
        aa_final_ter = os.path.join(output_dir, "backmapped_aa_final_ter.pdb")
        remove_ot2_atoms(aa_final_with_ter_and_ot2, aa_final_ter)
        final_output = aa_final_ter
    else:
        # Add TER records even if no chain IDs (save intermediate with OT2 atoms)
        aa_final_with_ter_and_ot2 = os.path.join(output_dir, "backmapped_aa_final_with_OT2.pdb")
        add_ter_records(aa_box_pdb, aa_final_with_ter_and_ot2)
        print(f"Saved intermediate file with OT2 atoms: {aa_final_with_ter_and_ot2}")
        
        # Remove OT2 atoms to create final output
        aa_final_ter = os.path.join(output_dir, "backmapped_aa_final_ter.pdb")
        remove_ot2_atoms(aa_final_with_ter_and_ot2, aa_final_ter)
        final_output = aa_final_ter

    print("="*60)
    print("Backmapping completed successfully!")
    print(f"Final structure (without OT2): {final_output}")
    print("="*60)

if __name__ == "__main__":
    run_backmapping()
