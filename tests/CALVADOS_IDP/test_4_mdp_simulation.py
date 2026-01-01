#!/usr/bin/env python3
"""
Test 4: MDP Simulation with config_mdp.yaml
使用 config_mdp.yaml 进行 MDP（基于结构的折叠蛋白）模拟测试。
"""

import os
import sys
from pathlib import Path


def test_config_loading():
    """测试 4.1: 配置文件加载"""
    print("\n" + "=" * 60)
    print("测试 4.1: MDP 配置加载 (config_mdp.yaml)")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulationConfig, CGSimulator, ComponentType, TopologyType
        
        config_path = Path(__file__).parent / "config_mdp.yaml"
        
        if not config_path.exists():
            print(f"✗ 配置文件不存在: {config_path}")
            return None
        
        config = CGSimulationConfig.from_yaml(str(config_path))
        
        print(f"✓ 配置文件加载成功")
        print(f"  - 系统: {config.system_name}")
        print(f"  - 盒子: {config.box}")
        print(f"  - 拓扑: {config.topol}")
        print(f"  - 温度: {config.temperature} K")
        print(f"  - 离子强度: {config.ionic} M")
        
        # 验证 CUBIC 拓扑
        if config.topol != TopologyType.CUBIC:
            print(f"✗ 期望 CUBIC 拓扑，实际 {config.topol}")
            return None
        print(f"  ✓ CUBIC 拓扑验证通过")
        
        return config
        
    except Exception as e:
        print(f"✗ 配置加载失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_config_validation(config):
    """测试 4.2: 配置验证"""
    print("\n" + "=" * 60)
    print("测试 4.2: MDP 配置验证")
    print("=" * 60)
    
    from multiscale2.src.cg import ComponentType
    
    try:
        errors = config.validate()
        if errors:
            print(f"✗ 配置验证失败: {errors}")
            return False
        print(f"✓ 配置验证通过")
        
        # 验证组件
        if len(config.components) != 1:
            print(f"✗ 期望 1 个组件，实际 {len(config.components)}")
            return False
        
        comp = config.components[0]
        print(f"  - 组件: {comp.name}")
        print(f"  - 类型: {comp.type.value}")
        print(f"  - 分子数: {comp.nmol}")
        print(f"  - PDB 文件: {comp.fpdb}")
        print(f"  - 约束: {comp.restraint}")
        
        if comp.type != ComponentType.MDP:
            print(f"✗ 期望 MDP 组件，实际 {comp.type}")
            return False
        if not comp.restraint:
            print(f"✗ MDP 应启用结构约束")
            return False
        if 'TDP43:' not in str(comp.fdomains):
            print(f"✗ fdomains 应包含 TDP43 定义")
            return False
        print(f"  ✓ MDP 组件验证通过")
        
        return True
        
    except Exception as e:
        print(f"✗ 配置验证失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_pdb_loading(config):
    """测试 4.3: PDB 文件加载并提取残基数"""
    print("\n" + "=" * 60)
    print("测试 4.3: PDB 文件加载")
    print("=" * 60)
    
    from multiscale2.src.cg import ComponentType
    
    try:
        from Bio import PDB
        
        for comp in config.components:
            if comp.type == ComponentType.MDP:
                if comp.fpdb and os.path.exists(comp.fpdb):
                    # 使用 Biopython 解析 PDB
                    parser = PDB.PDBParser(QUIET=True)
                    structure = parser.get_structure(comp.name, comp.fpdb)
                    
                    # 计算残基数
                    nres = 0
                    for model in structure:
                        for chain in model:
                            nres += len(list(chain))
                    
                    comp.nres = nres
                    print(f"  ✓ {comp.name}: {comp.fpdb}")
                    print(f"    - 残基数: {nres}")
                else:
                    print(f"  ⚠ {comp.name}: PDB 文件不存在")
        
        print(f"✓ PDB 加载完成")
        return True
        
    except Exception as e:
        print(f"✗ PDB 加载失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_mdp_simulation_run(config):
    """测试 4.4: 运行 MDP 模拟 (使用 CGSimulator)"""
    print("\n" + "=" * 60)
    print("测试 4.4: 运行 MDP 模拟 (CGSimulator)")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulator
        import shutil
        
        # 输出目录
        test_output_dir = Path(__file__).parent / "test4_output"
        test_output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"输出目录: {test_output_dir}")
        print(f"  系统: {config.system_name}")
        print(f"  拓扑: {config.topol.value}")
        print(f"  分子数: {config.total_molecules()}")
        print(f"  步数: {config.simulation.steps}")
        print(f"  约束: {config.components[0].restraint}")
        
        # 创建 Simulator 并运行
        sim = CGSimulator(config)
        
        # Setup
        sim.setup(str(test_output_dir), overwrite=True)
        
        # 运行模拟
        result = sim.run_calvados(gpu_id=0)
        
        if result.success:
            print(f"\n✓ MDP 模拟完成!")
            print(f"  输出目录: {result.output_dir}")
            
            if result.trajectory and os.path.exists(result.trajectory):
                size_mb = os.path.getsize(result.trajectory) / 1024 / 1024
                print(f"  轨迹: {os.path.basename(result.trajectory)} ({size_mb:.1f} MB)")
            
            if result.structure and os.path.exists(result.structure):
                print(f"  结构: {os.path.basename(result.structure)}")
                print(f"  ✓ 直接使用 CALVADOS 生成的 PDB")
            
            # 验证 final.pdb 的 chain 编号
            if result.structure:
                from multiscale2.src import extract_coordinates_from_pdb
                coords = extract_coordinates_from_pdb(result.structure)
                print(f"  ✓ Chain IDs: {list(coords.keys())}")
                print(f"  ✓ Total chains: {len(coords)}")
                for chain, crds in coords.items():
                    print(f"    {chain}: {len(crds)} residues")
            
            return True
        else:
            print(f"\n✗ 模拟失败!")
            for error in result.errors:
                print(f"  Error: {error}")
            return False
        
    except Exception as e:
        print(f"\n✗ 模拟运行失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_coordinate_consistency():
    """测试 4.5: 验证 final.pdb 与 raw PDB 的坐标一致性 (使用 MDAnalysis)
    
    说明：
    - 时间戳 PDB 文件将所有 100 个 TDP43 分子放在 Chain A 中（残基 1-4140）
    - final.pdb 将每个分子分配到独立的链（Chain 0-99，每个链 414 个残基）
    - 两者坐标应该完全相同，只是组织方式不同
    """
    print("\n" + "=" * 60)
    print("测试 4.5: 坐标一致性验证 (MDAnalysis)")
    print("=" * 60)
    
    try:
        import numpy as np
        
        # 尝试导入 MDAnalysis
        try:
            import MDAnalysis as mda
        except ImportError:
            print("  ⚠ MDAnalysis 未安装，跳过此测试")
            print("  💡 安装命令: pip install MDAnalysis")
            return True
        
        # 获取输出目录
        test_output_dir = Path(__file__).parent / "test4_output"

        # final.pdb 在 TDP43_CG 子目录中
        cg_output_dir = test_output_dir / "TDP43_CG"
        final_pdb = cg_output_dir / "final.pdb"
        if not final_pdb.exists():
            print(f"✗ final.pdb 不存在: {final_pdb}")
            return False

        raw_dir = cg_output_dir / "raw"
        
        if not raw_dir.exists():
            print(f"✗ raw 目录不存在: {raw_dir}")
            return False
        
        # 查找时间戳 PDB 文件（模拟结束时的最终结构）
        raw_pdb = None
        for f in raw_dir.glob("*.pdb"):
            # 跳过 checkpoint.pdb、top.pdb、system.pdb
            if f.name in ('checkpoint.pdb', 'top.pdb', 'system.pdb'):
                continue
            # 跳过非时间戳文件
            if not any(c in f.name for c in ['_', '20']):
                continue
            raw_pdb = f
            break
        
        if raw_pdb is None or not raw_pdb.exists():
            print(f"✗ 时间戳 PDB 文件不存在")
            return False
        
        print(f"  比较文件:")
        print(f"    - final.pdb: {final_pdb.name}")
        print(f"    - raw PDB: {raw_pdb.name}")
        
        # 使用 MDAnalysis 读取两个结构
        u_final = mda.Universe(str(final_pdb))
        u_raw = mda.Universe(str(raw_pdb))
        
        # 获取两个结构的 CA 原子坐标
        atoms_final = u_final.select_atoms("name CA")
        atoms_raw = u_raw.select_atoms("name CA")
        
        print(f"  final.pdb: {len(atoms_final)} CA 原子 (100 分子 x 414 残基)")
        print(f"  raw PDB: {len(atoms_raw)} CA 原子 (1 链 x 4140 残基)")
        
        # 检查原子数是否一致
        if len(atoms_final) != len(atoms_raw):
            print(f"✗ 原子数不匹配: {len(atoms_final)} vs {len(atoms_raw)}")
            return False
        
        # 提取坐标（按位置顺序）
        coords_final = atoms_final.positions  # Angstrom
        coords_raw = atoms_raw.positions      # Angstrom
        
        # 计算偏差
        diff = coords_final - coords_raw
        distances = np.linalg.norm(diff, axis=1)
        
        max_dist = np.max(distances)
        mean_dist = np.mean(distances)
        rmsd = np.sqrt(np.mean(distances ** 2))
        
        # 设置阈值
        tolerance = 0.1  # Angstrom - 由于是相同坐标，应该几乎为 0
        
        print(f"\n  坐标偏差统计:")
        print(f"    - 最大偏差: {max_dist:.6f} Å")
        print(f"    - 平均偏差: {mean_dist:.6f} Å")
        print(f"    - RMSD: {rmsd:.6f} Å")
        print(f"    - 阈值: {tolerance} Å")
        
        # 检查超出阈值的比例
        outliers = np.sum(distances > tolerance)
        if outliers > 0:
            outlier_pct = 100.0 * outliers / len(distances)
            print(f"  ⚠ {outliers} 个原子 ({outlier_pct:.2f}%) 偏差超出阈值")
            if outlier_pct < 1.0:
                # 打印偏差最大的几个原子
                outlier_indices = np.argsort(distances)[-5:][::-1]
                for idx in outlier_indices:
                    print(f"    Atom {idx}: 偏差 = {distances[idx]:.6f} Å")
        
        # 判断是否通过
        if max_dist <= tolerance:
            print(f"\n  ✓ 坐标一致性验证通过!")
            print(f"    final.pdb 与 {raw_pdb.name} 的坐标完全一致")
            return True
        else:
            print(f"\n  ✗ 坐标一致性验证失败!")
            print(f"    最大偏差 {max_dist:.6f} Å 超出阈值 {tolerance} Å")
            return False
        
    except Exception as e:
        print(f"\n✗ 坐标验证失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """主测试函数"""
    print("=" * 60)
    print("MDP 模拟测试 (config_mdp.yaml) - CGSimulator")
    print("=" * 60)
    
    results = []
    
    # 基础导入测试
    print("运行基础导入测试...")
    try:
        from multiscale2.src import CGSimulationConfig, CGSimulator, ComponentType, TopologyType
        print("✓ CG 模块可用")
    except ImportError as e:
        print(f"✗ 导入失败: {e}")
        return False
    
    # 测试 4.1: 配置加载
    config = test_config_loading()
    if config is None:
        results.append(("配置加载", False))
    else:
        results.append(("配置加载", True))
    
    # 测试 4.2: 配置验证
    if config is not None:
        results.append(("配置验证", test_config_validation(config)))
    else:
        results.append(("配置验证", False))
    
    # 测试 4.3: PDB 加载
    if config is not None:
        results.append(("PDB 加载", test_pdb_loading(config)))
    else:
        results.append(("PDB 加载", False))
    
    # 测试 4.4: 运行模拟
    if config is not None:
        results.append(("模拟运行", test_mdp_simulation_run(config)))
    else:
        results.append(("模拟运行", False))
    
    # 测试 4.5: 坐标一致性验证 (如果模拟已运行)
    test_output_dir = Path(__file__).parent / "test4_output"
    cg_output_dir = test_output_dir / "TDP43_CG"
    if config is not None and cg_output_dir.exists():
        results.append(("坐标一致性验证", test_coordinate_consistency()))
    else:
        results.append(("坐标一致性验证", None))  # 跳过
    
    # 总结
    print("\n" + "=" * 60)
    print("测试结果总结")
    print("=" * 60)
    for name, result in results:
        if result is None:
            status = "⊘"
        elif result is True:
            status = "✓"
        else:
            status = "✗"
        print(f"  {status} {name}")
    
    # 只检查明确失败的测试
    failed_results = [r for r in results if r[1] is False]
    all_passed = len(failed_results) == 0
    print()
    if all_passed:
        print("🎉 所有 MDP 测试通过!")
        return 0
    else:
        print(f"⚠️  {len(failed_results)} 个测试失败")
        return 1


if __name__ == "__main__":
    sys.exit(main())
