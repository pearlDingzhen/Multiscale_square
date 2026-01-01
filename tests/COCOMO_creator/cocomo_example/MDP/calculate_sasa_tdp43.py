#!/usr/bin/env python3
"""
Calculate SASA for TDP43 and prepare surface file for simulation.
Uses mdsim library same as cocomo-main.
"""

import numpy as np
import os
from pathlib import Path

from mdsim import PDBReader


def calculate_sasa_tdp43(pdb_file, n_sphere_points=1920):
    """Calculate SASA for each residue using mdsim."""
    reader = PDBReader(pdb_file)
    model = reader[0]
    sasa = model.sasa_by_residue(n_sphere_points=n_sphere_points)
    return sasa


def prepare_surface_file(sasa_values, domains, output_file):
    """
    Prepare surface file for TDP43 simulation.
    
    Args:
        sasa_values: List of SASA values per residue (nm²)
        domains: List of [start, end] for each domain (1-based, inclusive)
        output_file: Path to output surface file
    """
    n_residues = len(sasa_values)
    surface_vector = np.array(sasa_values)
    
    # 根据config_mdp.yaml中的domain定义
    # Domain区域: 代码会将其设为5.0，所以我们只需要提供linker的值
    # Linker区域: 保持原始SASA值
    
    # 创建一个mask来标识哪些是linker区域
    is_linker = np.ones(n_residues, dtype=bool)
    
    for start, end in domains:
        # 转换为0-based索引
        start_0 = start - 1
        end_0 = end
        if start_0 >= 0 and end_0 <= n_residues:
            is_linker[start_0:end_0] = False
    
    # 对于linker区域，我们使用计算得到的SASA值
    # 对于domain区域，虽然代码会覆盖为5.0，但我们也写入原始值作为参考
    
    # 写入surface文件
    with open(output_file, 'w') as f:
        for i in range(n_residues):
            if is_linker[i]:
                # Linker区域: 写入计算得到的SASA
                f.write(f"{surface_vector[i]:.3f}\n")
            else:
                # Domain区域: 写入5.0（会被代码使用）
                f.write("5.000\n")
    
    print(f"Surface file written to: {output_file}")
    
    # 返回统计信息
    n_linker = np.sum(is_linker)
    n_domain = n_residues - n_linker
    
    return {
        'total_residues': n_residues,
        'n_linker': n_linker,
        'n_domain': n_domain,
        'linker_sasa_mean': np.mean(surface_vector[is_linker]) if n_linker > 0 else 0,
        'linker_sasa_std': np.std(surface_vector[is_linker]) if n_linker > 0 else 0,
    }


def main():
    # TDP43 domains from config_mdp.yaml (1-based, inclusive)
    domains = [
        [3, 76],
        [106, 176],
        [192, 260],
        [320, 334],
    ]
    
    pdb_file = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/TDP43.pdb'
    output_file = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/surface'
    
    print("="*70)
    print("TDP43 SASA Calculation and Surface File Preparation")
    print("="*70)
    
    # 计算SASA
    print(f"\n1. Loading structure from: {pdb_file}")
    print("2. Calculating SASA using mdsim...")
    sasa_values = calculate_sasa_tdp43(pdb_file, n_sphere_points=1920)
    
    print(f"\n3. TDP43 Domain Definition (from config_mdp.yaml):")
    for i, (start, end) in enumerate(domains, 1):
        length = end - start + 1
        print(f"   Domain {i}: residues {start} - {end} ({length} residues)")
    
    # 准备surface文件
    print(f"\n4. Preparing surface file...")
    stats = prepare_surface_file(sasa_values, domains, output_file)
    
    print(f"\n5. Summary:")
    print(f"   Total residues: {stats['total_residues']}")
    print(f"   Domain residues: {stats['n_domain']}")
    print(f"   Linker residues: {stats['n_linker']}")
    print(f"   Linker SASA mean: {stats['linker_sasa_mean']:.3f} nm²")
    print(f"   Linker SASA std: {stats['linker_sasa_std']:.3f} nm²")
    
    # 显示前10个和后10个值
    print(f"\n6. Surface file preview (first 10 and last 10 values):")
    with open(output_file, 'r') as f:
        lines = f.readlines()
        print("   First 10:")
        for i, line in enumerate(lines[:10], 1):
            print(f"     {i}: {line.strip()}")
        print("   ...")
        print("   Last 10:")
        for i, line in enumerate(lines[-10:], len(lines)-9):
            print(f"     {i}: {line.strip()}")
    
    print(f"\n" + "="*70)
    print("Surface file ready for simulation!")
    print("="*70)


if __name__ == '__main__':
    main()
