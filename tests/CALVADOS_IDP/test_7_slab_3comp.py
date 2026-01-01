#!/usr/bin/env python3
"""
Test 7: Three-Component SLAB Simulation (TDP43_CTD + FUS + TDP43)
使用 config_idp_slab_3comp.yaml 进行三组分 SLAB 相分离模拟测试。
"""

import os
import sys
from pathlib import Path


def test_config_loading():
    """测试 7.1: 配置文件加载"""
    print("\n" + "=" * 60)
    print("测试 7.1: 三组分 SLAB 配置加载")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulationConfig, CGSimulator, ComponentType
        
        config_path = Path(__file__).parent / "config_idp_slab_3comp.yaml"
        
        if not config_path.exists():
            print(f"✗ 配置文件不存在: {config_path}")
            return None
        
        config = CGSimulationConfig.from_yaml(str(config_path))
        
        print(f"✓ 配置文件加载成功")
        print(f"  - 系统: {config.system_name}")
        print(f"  - 拓扑: {config.topol}")
        print(f"  - 组件数: {len(config.components)}")
        
        # 验证组件
        if len(config.components) != 3:
            print(f"✗ 期望 3 个组件，实际 {len(config.components)}")
            return None
        
        if config.components[0].name != "TDP43_CTD":
            print(f"✗ 期望 TDP43_CTD，实际 {config.components[0].name}")
            return None
        if config.components[1].name != "FUS":
            print(f"✗ 期望 FUS，实际 {config.components[1].name}")
            return None
        if config.components[2].name != "TDP43":
            print(f"✗ 期望 TDP43，实际 {config.components[2].name}")
            return None
        print(f"  ✓ 组件验证通过")
        
        # 验证组件类型
        if config.components[0].type != ComponentType.IDP:
            print(f"✗ 第一个组件应为 IDP，实际 {config.components[0].type}")
            return None
        if config.components[1].type != ComponentType.IDP:
            print(f"✗ 第二个组件应为 IDP，实际 {config.components[1].type}")
            return None
        if config.components[2].type != ComponentType.MDP:
            print(f"✗ 第三个组件应为 MDP，实际 {config.components[2].type}")
            return None
        print(f"  ✓ 组件类型验证通过")
        
        return config
        
    except Exception as e:
        print(f"✗ 配置加载失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_config_validation(config):
    """测试 7.2: 配置验证"""
    print("\n" + "=" * 60)
    print("测试 7.2: 三组分 SLAB 配置验证")
    print("=" * 60)
    
    from multiscale2.src.cg import ComponentType
    
    try:
        errors = config.validate()
        if errors:
            print(f"✗ 配置验证失败: {errors}")
            return False
        print(f"✓ 配置验证通过")
        total_mol = sum(c.nmol for c in config.components)
        print(f"  - 总分子数: {total_mol}")
        return True
        
    except Exception as e:
        print(f"✗ 配置验证失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_fasta_and_pdb_loading(config):
    """测试 7.3: FASTA 和 PDB 文件加载"""
    print("\n" + "=" * 60)
    print("测试 7.3: FASTA 和 PDB 文件加载")
    print("=" * 60)
    
    try:
        from Bio import SeqIO
        from Bio import PDB
        from multiscale2.src.cg import ComponentType
        
        # 处理 IDP 组件 - FASTA
        fasta_path = Path(__file__).parent / "TDP43_CTD.fasta"
        if fasta_path.exists():
            records = list(SeqIO.parse(str(fasta_path), "fasta"))
            fasta_dict = {rec.id: str(rec.seq) for rec in records}
            
            for comp in config.components:
                if comp.type == ComponentType.IDP:
                    if comp.name in fasta_dict:
                        comp.seq = fasta_dict[comp.name]
                        print(f"  ✓ {comp.name}: {len(comp.seq)} residues")
                    else:
                        comp.seq = fasta_dict.get("TDP43_CTD", "")
                        if comp.seq:
                            print(f"  ⚠ {comp.name}: 使用 TDP43_CTD 序列 ({len(comp.seq)} residues)")
        
        # 处理 MDP 组件 - PDB
        for comp in config.components:
            if comp.type == ComponentType.MDP:
                if comp.fpdb and os.path.exists(comp.fpdb):
                    parser = PDB.PDBParser(QUIET=True)
                    structure = parser.get_structure(comp.name, comp.fpdb)
                    
                    nres = 0
                    for model in structure:
                        for chain in model:
                            nres += len(list(chain))
                    
                    comp.nres = nres
                    print(f"  ✓ {comp.name}: {comp.fpdb}")
                    print(f"    - 残基数: {nres}")
                else:
                    print(f"  ⚠ {comp.name}: PDB 文件不存在")
        
        print(f"✓ 文件加载完成")
        return True
        
    except Exception as e:
        print(f"✗ 文件加载失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_slab_simulation_run(config):
    """测试 7.4: 运行三组分 SLAB 模拟 (使用 CGSimulator)"""
    print("\n" + "=" * 60)
    print("测试 7.4: 运行 SLAB 三组分模拟 (CGSimulator)")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulator
        import shutil
        
        # 输出目录
        test_output_dir = Path(__file__).parent / "test7_output"
        test_output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"输出目录: {test_output_dir}")
        print(f"  系统: {config.system_name}")
        print(f"  拓扑: {config.topol.value}")
        print(f"  分子数: {config.total_molecules()}")
        print(f"  步数: {config.simulation.steps}")
        
        # 创建 Simulator 并运行
        sim = CGSimulator(config)
        
        # Setup
        sim.setup(str(test_output_dir), overwrite=True)
        
        # 运行模拟
        result = sim.run_calvados(gpu_id=0)
        
        if result.success:
            print(f"\n✓ SLAB 三组分模拟完成!")
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


def main():
    """主测试函数"""
    print("=" * 60)
    print("三组分 SLAB 模拟测试 (TDP43_CTD + FUS + TDP43) - CGSimulator")
    print("=" * 60)
    
    results = []
    
    # 基础导入测试
    print("运行基础导入测试...")
    try:
        from multiscale2.src import CGSimulationConfig, CGSimulator, ComponentType
        print("✓ CG 模块可用")
    except ImportError as e:
        print(f"✗ 导入失败: {e}")
        return False
    
    # 测试 7.1: 配置加载
    config = test_config_loading()
    if config is None:
        results.append(("配置加载", False))
    else:
        results.append(("配置加载", True))
    
    # 测试 7.2: 配置验证
    if config is not None:
        results.append(("配置验证", test_config_validation(config)))
    else:
        results.append(("配置验证", False))
    
    # 测试 7.3: 文件加载
    if config is not None:
        results.append(("文件加载", test_fasta_and_pdb_loading(config)))
    else:
        results.append(("文件加载", False))
    
    # 测试 7.4: 运行模拟
    if config is not None:
        results.append(("模拟运行", test_slab_simulation_run(config)))
    else:
        results.append(("模拟运行", False))
    
    # 总结
    print("\n" + "=" * 60)
    print("测试结果总结")
    print("=" * 60)
    for name, result in results:
        status = "✓" if result is True else "✗"
        print(f"  {status} {name}")
    
    all_passed = all(r is True for _, r in results)
    print()
    if all_passed:
        print("🎉 所有三组分 SLAB 测试通过!")
        return 0
    else:
        print("⚠️  部分测试失败")
        return 1


if __name__ == "__main__":
    sys.exit(main())
