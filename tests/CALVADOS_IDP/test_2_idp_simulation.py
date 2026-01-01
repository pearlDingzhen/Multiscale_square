#!/usr/bin/env python3
"""
Test 2: IDP Simulation with config_idp.yaml
使用 config_idp.yaml 进行 IDP 模拟测试

从配置文件加载配置，运行 CALVADOS 模拟。
"""

import os
import sys
from pathlib import Path


def test_config_validation():
    """测试配置文件验证"""
    print("=" * 60)
    print("测试 2.1: 配置文件验证")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulationConfig
        
        config_path = Path(__file__).parent / "config_idp.yaml"
        
        if not config_path.exists():
            print(f"✗ 配置文件不存在: {config_path}")
            return False
        
        config = CGSimulationConfig.from_yaml(str(config_path))
        
        print(f"✓ 配置文件加载成功")
        print(f"  - 系统: {config.system_name}")
        print(f"  - 盒子: {config.box}")
        print(f"  - 拓扑: {config.topol.value}")
        
        # 验证
        errors = config.validate()
        if errors:
            print(f"✗ 验证失败: {errors}")
            return False
        
        print(f"  ✓ 配置验证通过")
        return True
        
    except Exception as e:
        print(f"✗ 配置验证失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_wrapper_creation():
    """测试 Wrapper 创建"""
    print("\n" + "=" * 60)
    print("测试 2.2: CALVADOS Wrapper 创建")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulationConfig
        from multiscale2.src.calvados_wrapper import CalvadosWrapper
        
        config_path = Path(__file__).parent / "config_idp.yaml"
        config = CGSimulationConfig.from_yaml(str(config_path))
        
        wrapper = CalvadosWrapper(config)
        
        print(f"✓ CALVADOS Wrapper 创建成功")
        return True
        
    except Exception as e:
        print(f"✗ Wrapper 创建失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_config_file_generation():
    """测试配置文件生成"""
    print("\n" + "=" * 60)
    print("测试 2.3: 配置文件生成")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulationConfig
        from multiscale2.src.calvados_wrapper import CalvadosWrapper
        import tempfile
        import shutil
        
        config_path = Path(__file__).parent / "config_idp.yaml"
        config = CGSimulationConfig.from_yaml(str(config_path))
        
        wrapper = CalvadosWrapper(config)
        
        # 使用临时目录
        output_dir = Path(tempfile.mkdtemp(prefix="idp_test_"))
        
        try:
            files = wrapper.write(str(output_dir), overwrite=True)
            
            print(f"✓ 配置文件生成成功")
            print(f"  - config.yaml: {files['config']}")
            print(f"  - components.yaml: {files['components']}")
            
            # 验证文件存在
            if os.path.exists(files['config']) and os.path.exists(files['components']):
                print(f"  ✓ 文件验证通过")
                return True
            else:
                print(f"✗ 文件不存在")
                return False
            
        finally:
            # 清理临时目录
            if output_dir.exists():
                shutil.rmtree(output_dir)
        
    except Exception as e:
        print(f"✗ 配置文件生成失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_simulation_run():
    """运行模拟测试"""
    print("\n" + "=" * 60)
    print("测试 2.4: 运行 CALVADOS 模拟 (CGSimulator)")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulationConfig, CGSimulator
        from Bio import SeqIO
        import shutil
        
        config_path = Path(__file__).parent / "config_idp.yaml"
        config = CGSimulationConfig.from_yaml(str(config_path))
        
        # 解析 FASTA 并填充 comp.seq
        fasta_path = Path(__file__).parent / "TDP43_CTD.fasta"
        if fasta_path.exists():
            records = list(SeqIO.parse(str(fasta_path), "fasta"))
            fasta_dict = {rec.id: str(rec.seq) for rec in records}
            for comp in config.components:
                if comp.type.value == "idp":
                    if comp.name in fasta_dict:
                        comp.seq = fasta_dict[comp.name]
                        print(f"  序列: {comp.name} = {len(comp.seq)} residues")
        
        # 输出目录 - 使用 overwrite=True，无需手动删除
        test_output_dir = Path(__file__).parent / "test2_output"
        test_output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"输出目录: {test_output_dir}")
        print(f"  系统: {config.system_name}")
        print(f"  分子数: {config.total_molecules()}")
        print(f"  步数: {config.simulation.steps}")
        
        # 创建 Simulator 并运行
        sim = CGSimulator(config)
        
        # Setup
        sim.setup(str(test_output_dir), overwrite=True)
        
        # 运行模拟
        result = sim.run_calvados(gpu_id=0)
        
        if result.success:
            print(f"\n✓ 模拟完成!")
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
    print("\n" + "=" * 60)
    print("IDP 模拟测试 (config_idp.yaml)")
    print("=" * 60)
    print()
    
    results = []
    
    # 先运行导入测试
    print("运行基础导入测试...")
    from multiscale2.src import CGSimulationConfig
    print("✓ CG 模块可用")
    print()
    
    results.append(("配置验证", test_config_validation()))
    results.append(("Wrapper 创建", test_wrapper_creation()))
    results.append(("文件生成", test_config_file_generation()))
    results.append(("模拟运行", test_simulation_run()))
    
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
        print("🎉 所有 IDP 测试通过!")
        return 0
    else:
        print("⚠️  部分测试失败")
        return 1


if __name__ == "__main__":
    sys.exit(main())

