#!/usr/bin/env python3
"""
Test 7: 三组分 SLAB 模拟 (TDP43_CTD + FUS + TDP43)

测试 IDP + IDP + MDP 三组分混合模拟在 SLAB 拓扑下的运行。
"""

import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from multiscale2.src.cg import CGSimulationConfig, CGComponent, ComponentType, TopologyType
from multiscale2.src.calvados_wrapper import CalvadosWrapper

def test_config_loading():
    """测试 7.1: 配置文件加载"""
    print("\n============================================================")
    print("测试 7.1: 三组分 SLAB 配置加载")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab_3comp.yaml"
    assert config_path.exists(), f"配置文件不存在: {config_path}"
    
    config = CGSimulationConfig.from_yaml(str(config_path))
    print(f"✓ 配置文件加载成功")
    print(f"  - 系统: {config.system_name}")
    print(f"  - 拓扑: {config.topol}")
    print(f"  - 组件数: {len(config.components)}")
    
    # 验证组件
    assert len(config.components) == 3, f"期望 3 个组件，实际 {len(config.components)}"
    assert config.components[0].name == "TDP43_CTD"
    assert config.components[1].name == "FUS"
    assert config.components[2].name == "TDP43"
    print(f"  ✓ 组件验证通过")
    
    # 验证组件类型
    assert config.components[0].type == ComponentType.IDP
    assert config.components[1].type == ComponentType.IDP
    assert config.components[2].type == ComponentType.MDP
    print(f"  ✓ 组件类型验证通过")
    
    return config

def test_config_validation(config):
    """测试 7.2: 配置验证"""
    print("\n============================================================")
    print("测试 7.2: 三组分 SLAB 配置验证")
    print("============================================================")
    
    errors = config.validate()
    assert len(errors) == 0, f"配置验证失败: {errors}"
    print(f"✓ 配置验证通过")
    total_mol = sum(c.nmol for c in config.components)
    print(f"  - 总分子数: {total_mol}")

def test_slab_config_generation():
    """测试 7.3: SLAB 配置文件生成"""
    print("\n============================================================")
    print("测试 7.3: 三组分 SLAB 配置文件生成")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab_3comp.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))
    
    wrapper = CalvadosWrapper(config)
    
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = os.path.join(tmpdir, config.system_name + "_CG")
        wrapper.write(output_dir, overwrite=True)
        
        config_file = os.path.join(output_dir, "config.yaml")
        components_file = os.path.join(output_dir, "components.yaml")
        
        assert os.path.exists(config_file), f"config.yaml 不存在"
        assert os.path.exists(components_file), f"components.yaml 不存在"
        
        print(f"✓ 配置文件生成成功")
        
        import yaml
        with open(config_file) as f:
            calvados_config = yaml.safe_load(f)
        
        assert calvados_config.get('slab_width') is not None
        assert calvados_config.get('topol') == 'slab'
        assert calvados_config['slab_width'] == 40.0
        
        print(f"  ✓ SLAB 参数验证通过")
        
        # 验证所有组件
        with open(components_file) as f:
            calvados_components = yaml.safe_load(f)
        
        system = calvados_components.get('system', {})
        assert 'TDP43_CTD' in system
        assert 'FUS' in system
        assert 'TDP43' in system
        print(f"  ✓ 所有三组件验证通过")

def test_slab_simulation_run():
    """测试 7.5: 运行 SLAB 三组分模拟"""
    print("\n============================================================")
    print("测试 7.5: 运行 SLAB 三组分模拟 (TDP43_CTD + FUS + TDP43)")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab_3comp.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))
    
    output_base = Path(__file__).parent
    output_dir = str(output_base)
    
    # 清理旧目录
    import shutil
    cg_output_dir = os.path.join(output_dir, f"{config.system_name}_CG")
    if os.path.exists(cg_output_dir):
        shutil.rmtree(cg_output_dir)
    
    wrapper = CalvadosWrapper(config)
    wrapper.run(output_dir)
    
    # 验证输出
    assert os.path.exists(cg_output_dir), f"输出目录不存在: {cg_output_dir}"
    assert os.path.exists(os.path.join(cg_output_dir, "final.pdb")), "缺少 final.pdb"
    assert os.path.exists(os.path.join(cg_output_dir, "trajectory.dcd")), "缺少 trajectory.dcd"
    
    raw_dir = os.path.join(cg_output_dir, "raw")
    assert os.path.exists(raw_dir), f"raw 目录不存在"
    
    print(f"✓ SLAB 三组分模拟完成!")
    print(f"  输出目录: {cg_output_dir}")
    print(f"  - final.pdb: {os.path.exists(os.path.join(cg_output_dir, 'final.pdb'))}")
    print(f"  - trajectory.dcd: {os.path.exists(os.path.join(cg_output_dir, 'trajectory.dcd'))}")

def main():
    """主测试函数"""
    print("=" * 60)
    print("三组分 SLAB 模拟测试 (TDP43_CTD + FUS + TDP43)")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulationConfig
        from multiscale2.src.calvados_wrapper import CalvadosWrapper
        print("✓ CG 模块可用")
    except ImportError as e:
        print(f"✗ 导入失败: {e}")
        return False
    
    try:
        config = test_config_loading()
        test_config_validation(config)
        test_slab_config_generation()
        test_slab_simulation_run()
        
        print("\n" + "=" * 60)
        print("🎉 所有三组分 SLAB 测试通过!")
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
