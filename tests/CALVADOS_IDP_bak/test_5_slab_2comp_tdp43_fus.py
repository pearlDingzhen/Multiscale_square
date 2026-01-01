#!/usr/bin/env python3
"""
Test 5: 双组分 SLAB 模拟 (TDP43_CTD + FUS)

测试多组分 IDP + IDP 模拟在 SLAB 拓扑下的运行。
"""

import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from multiscale2.src.cg import CGSimulationConfig, CGComponent, ComponentType, TopologyType
from multiscale2.src.calvados_wrapper import CalvadosWrapper

def test_config_loading():
    """测试 5.1: 配置文件加载"""
    print("\n============================================================")
    print("测试 5.1: 双组分 SLAB 配置加载 (TDP43_CTD + FUS)")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab_2comp_tdp43_fus.yaml"
    assert config_path.exists(), f"配置文件不存在: {config_path}"
    
    config = CGSimulationConfig.from_yaml(str(config_path))
    print(f"✓ 配置文件加载成功")
    print(f"  - 系统: {config.system_name}")
    print(f"  - 拓扑: {config.topol}")
    print(f"  - 组件数: {len(config.components)}")
    
    # 验证组件
    assert len(config.components) == 2, f"期望 2 个组件，实际 {len(config.components)}"
    assert config.components[0].name == "TDP43_CTD"
    assert config.components[1].name == "SHORT"
    print(f"  ✓ 组件验证通过")
    
    return config

def test_config_validation(config):
    """测试 5.2: 配置验证"""
    print("\n============================================================")
    print("测试 5.2: 双组分 SLAB 配置验证")
    print("============================================================")
    
    errors = config.validate()
    assert len(errors) == 0, f"配置验证失败: {errors}"
    print(f"✓ 配置验证通过")
    print(f"  - TDP43_CTD: type={config.components[0].type.value}, nmol={config.components[0].nmol}")
    print(f"  - FUS: type={config.components[1].type.value}, nmol={config.components[1].nmol}")

def test_slab_config_generation():
    """测试 5.3: SLAB 配置文件生成"""
    print("\n============================================================")
    print("测试 5.3: 双组分 SLAB 配置文件生成")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab_2comp_tdp43_fus.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))
    
    # 创建 wrapper
    wrapper = CalvadosWrapper(config)
    
    # 生成配置文件到临时目录
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
        print(f"  - config.yaml: {config_file}")
        print(f"  - components.yaml: {components_file}")
        
        # 读取并验证 SLAB 设置
        import yaml
        with open(config_file) as f:
            calvados_config = yaml.safe_load(f)
        
        # SLAB 参数在顶级键中
        assert calvados_config.get('slab_width') is not None, "缺少 slab_width"
        assert calvados_config.get('slab_eq') is not None, "缺少 slab_eq"
        assert calvados_config.get('topol') == 'slab', "topol 应为 slab"
        
        # 验证 SLAB 自动设置
        assert calvados_config['slab_eq'] == False, "slab_eq 应为 false"
        assert calvados_config['slab_width'] == 40.0, f"slab_width 应为 40.0 (box_z/2)，实际 {calvados_config['slab_width']}"
        
        print(f"  ✓ SLAB 参数验证通过")
        print(f"    - topol: {calvados_config['topol']}")
        print(f"    - slab_eq: {calvados_config['slab_eq']}")
        print(f"    - slab_width: {calvados_config['slab_width']}")

def test_slab_simulation_run():
    """测试 5.4: 运行 SLAB 双组分模拟"""
    print("\n============================================================")
    print("测试 5.4: 运行 SLAB 双组分模拟 (TDP43_CTD + FUS)")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab_2comp_tdp43_fus.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))
    
    # 输出到测试目录（wrapper 会自动添加 _CG 后缀）
    output_base = Path(__file__).parent
    output_dir = str(output_base)
    
    # 清理旧目录
    import shutil
    cg_output_dir = os.path.join(output_dir, f"{config.system_name}_CG")
    if os.path.exists(cg_output_dir):
        shutil.rmtree(cg_output_dir)
    
    # 创建 wrapper 并运行
    wrapper = CalvadosWrapper(config)
    wrapper.run(output_dir)
    
    # 验证输出
    assert os.path.exists(cg_output_dir), f"输出目录不存在: {cg_output_dir}"
    
    # 检查统一输出文件
    assert os.path.exists(os.path.join(cg_output_dir, "final.pdb")), "缺少 final.pdb"
    assert os.path.exists(os.path.join(cg_output_dir, "trajectory.dcd")), "缺少 trajectory.dcd"
    assert os.path.exists(os.path.join(cg_output_dir, "simulation.log")), "缺少 simulation.log"
    
    # 检查 raw 目录
    raw_dir = os.path.join(cg_output_dir, "raw")
    assert os.path.exists(raw_dir), f"raw 目录不存在: {raw_dir}"
    
    print(f"✓ SLAB 双组分模拟完成!")
    print(f"  输出目录: {cg_output_dir}")
    print(f"  - final.pdb: {os.path.exists(os.path.join(cg_output_dir, 'final.pdb'))}")
    print(f"  - trajectory.dcd: {os.path.exists(os.path.join(cg_output_dir, 'trajectory.dcd'))}")
    print(f"  - raw/: {os.listdir(raw_dir)}")

def main():
    """主测试函数"""
    print("=" * 60)
    print("双组分 SLAB 模拟测试 (TDP43_CTD + FUS)")
    print("=" * 60)
    
    # 基础导入测试
    try:
        from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType
        from multiscale2.src.calvados_wrapper import CalvadosWrapper
        print("✓ CG 模块可用")
    except ImportError as e:
        print(f"✗ 导入失败: {e}")
        return False
    
    try:
        # 测试 5.1: 配置加载
        config = test_config_loading()
        
        # 测试 5.2: 配置验证
        test_config_validation(config)
        
        # 测试 5.3: 配置文件生成
        test_slab_config_generation()
        
        # 测试 5.4: 运行模拟
        test_slab_simulation_run()
        
        print("\n" + "=" * 60)
        print("🎉 所有双组分 SLAB 测试通过!")
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
