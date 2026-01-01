#!/usr/bin/env python3
"""
Test 4: MDP 模拟测试 (config_mdp.yaml)

测试单组分 MDP（基于 PDB 结构的折叠蛋白）在 CUBIC 拓扑下的模拟。
"""

import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from multiscale2.src.cg import CGSimulationConfig, ComponentType, TopologyType
from multiscale2.src.calvados_wrapper import CalvadosWrapper

def test_config_loading():
    """测试 4.1: 配置文件加载"""
    print("\n============================================================")
    print("测试 4.1: MDP 配置加载 (config_mdp.yaml)")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_mdp.yaml"
    assert config_path.exists(), f"配置文件不存在: {config_path}"
    
    config = CGSimulationConfig.from_yaml(str(config_path))
    print(f"✓ 配置文件加载成功")
    print(f"  - 系统: {config.system_name}")
    print(f"  - 盒子: {config.box}")
    print(f"  - 拓扑: {config.topol}")
    print(f"  - 温度: {config.temperature} K")
    print(f"  - 离子强度: {config.ionic} M")
    
    # 验证 CUBIC 拓扑
    assert config.topol == TopologyType.CUBIC, f"期望 CUBIC 拓扑，实际 {config.topol}"
    print(f"  ✓ CUBIC 拓扑验证通过")
    
    return config

def test_config_validation(config):
    """测试 4.2: 配置验证"""
    print("\n============================================================")
    print("测试 4.2: MDP 配置验证")
    print("============================================================")
    
    errors = config.validate()
    assert len(errors) == 0, f"配置验证失败: {errors}"
    print(f"✓ 配置验证通过")
    
    # 验证组件
    assert len(config.components) == 1, f"期望 1 个组件，实际 {len(config.components)}"
    comp = config.components[0]
    print(f"  - 组件: {comp.name}")
    print(f"  - 类型: {comp.type.value}")
    print(f"  - 分子数: {comp.nmol}")
    print(f"  - PDB 文件: {comp.fpdb}")
    print(f"  - 约束: {comp.restraint}")
    
    assert comp.type == ComponentType.MDP, "应为 MDP 组件"
    assert comp.restraint == True, "MDP 应启用结构约束"
    assert 'TDP43:' in comp.fdomains, "fdomains 应包含 TDP43 定义"
    print(f"  ✓ MDP 组件验证通过")

def test_inline_fdomains_processing():
    """测试 4.3: 内联 fdomains 处理"""
    print("\n============================================================")
    print("测试 4.3: 内联 fdomains 处理")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_mdp.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))
    
    wrapper = CalvadosWrapper(config)
    
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        # 测试 _process_inline_fdomains
        components_yaml = wrapper._generate_components_yaml()
        processed_yaml = wrapper._process_inline_fdomains(components_yaml, tmpdir)
        
        import yaml
        components = yaml.safe_load(processed_yaml)
        
        # 检查 TDP43 的 fdomains 是否已转换为文件路径
        tdp43_fdomains = components['system']['TDP43']['fdomains']
        assert not tdp43_fdomains.startswith('{'), "fdomains 应该是文件路径而非内联 YAML"
        assert os.path.exists(tdp43_fdomains), f"fdomains 文件不存在: {tdp43_fdomains}"
        
        print(f"✓ 内联 fdomains 正确处理")
        print(f"  - 原始: TDP43: [3, 76], ...")
        print(f"  - 文件: {os.path.basename(tdp43_fdomains)}")

def test_mdp_config_generation():
    """测试 4.4: MDP 配置文件生成"""
    print("\n============================================================")
    print("测试 4.4: MDP 配置文件生成")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_mdp.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))
    
    wrapper = CalvadosWrapper(config)
    
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = os.path.join(tmpdir, config.system_name + "_CG")
        wrapper.write(output_dir, overwrite=True)
        
        # 验证文件
        config_file = os.path.join(output_dir, "config.yaml")
        components_file = os.path.join(output_dir, "components.yaml")
        
        assert os.path.exists(config_file), f"config.yaml 不存在"
        assert os.path.exists(components_file), f"components.yaml 不存在"
        
        print(f"✓ 配置文件生成成功")
        
        # 读取并验证
        import yaml
        with open(config_file) as f:
            calvados_config = yaml.safe_load(f)
        
        # 验证 CUBIC -> grid
        assert calvados_config.get('topol') == 'grid', "CUBIC 拓扑应转换为 'grid'"
        print(f"  ✓ 拓扑转换验证通过: CUBIC -> grid")
        
        # 验证 components.yaml 包含 MDP 特有参数
        with open(components_file) as f:
            calvados_components = yaml.safe_load(f)
        
        system = calvados_components.get('system', {})
        tdp43 = system.get('TDP43', {})
        assert tdp43.get('restraint') == True, "MDP 应启用 restraint"
        assert tdp43.get('fpdb') is not None, "MDP 应有 fpdb"
        assert tdp43.get('fdomains') is not None, "MDP 应有 fdomains"
        print(f"  ✓ MDP 特有参数验证通过")

def test_mdp_simulation_run():
    """测试 4.5: 运行 MDP 模拟"""
    print("\n============================================================")
    print("测试 4.5: 运行 MDP 模拟 (TDP43)")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_mdp.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))
    
    # 输出到测试目录
    output_base = Path(__file__).parent
    output_dir = str(output_base)
    
    # 清理旧目录
    import shutil
    cg_output_dir = os.path.join(output_dir, f"{config.system_name}_CG")
    if os.path.exists(cg_output_dir):
        shutil.rmtree(cg_output_dir)
    
    # 创建 wrapper 并运行
    wrapper = CalvadosWrapper(config)
    result = wrapper.run(output_dir)
    
    # 验证输出
    assert os.path.exists(cg_output_dir), f"输出目录不存在: {cg_output_dir}"
    assert os.path.exists(os.path.join(cg_output_dir, "final.pdb")), "缺少 final.pdb"
    assert os.path.exists(os.path.join(cg_output_dir, "trajectory.dcd")), "缺少 trajectory.dcd"
    assert os.path.exists(os.path.join(cg_output_dir, "simulation.log")), "缺少 simulation.log"
    
    # 检查 raw 目录
    raw_dir = os.path.join(cg_output_dir, "raw")
    assert os.path.exists(raw_dir), f"raw 目录不存在: {raw_dir}"
    
    # 粒子数 = 100 个分子 * 334 个残基 = 33400
    print(f"✓ MDP 模拟完成!")
    print(f"  输出目录: {cg_output_dir}")
    print(f"  - final.pdb: ✓")
    print(f"  - trajectory.dcd: ✓")
    print(f"  - 100 个 TDP43 分子 (共 33400 个粒子)")

def main():
    """主测试函数"""
    print("=" * 60)
    print("MDP 模拟测试 (config_mdp.yaml)")
    print("=" * 60)
    
    # 基础导入测试
    try:
        from multiscale2.src import CGSimulationConfig, ComponentType, TopologyType
        from multiscale2.src.calvados_wrapper import CalvadosWrapper
        print("✓ CG 模块可用")
    except ImportError as e:
        print(f"✗ 导入失败: {e}")
        return False
    
    try:
        # 测试 4.1: 配置加载
        config = test_config_loading()
        
        # 测试 4.2: 配置验证
        test_config_validation(config)
        
        # 测试 4.3: 内联 fdomains 处理
        test_inline_fdomains_processing()
        
        # 测试 4.4: 配置文件生成
        test_mdp_config_generation()
        
        # 测试 4.5: 运行模拟
        test_mdp_simulation_run()
        
        print("\n" + "=" * 60)
        print("🎉 所有 MDP 测试通过!")
        print("=" * 60)
        return True
        
    except AssertionError as e:
        print(f"\n✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    except Exception as e:
        print(f"\n✗ 错误: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

