#!/usr/bin/env python3
"""
Test 3: SLAB 模拟测试 (config_idp_slab.yaml)

测试单组分 IDP 在 SLAB 拓扑下的相分离模拟。
"""

import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from multiscale2.src.cg import CGSimulationConfig, ComponentType, TopologyType
from multiscale2.src.calvados_wrapper import CalvadosWrapper

def test_config_loading():
    """测试 3.1: 配置文件加载"""
    print("\n============================================================")
    print("测试 3.1: SLAB 配置加载 (config_idp_slab.yaml)")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab.yaml"
    assert config_path.exists(), f"配置文件不存在: {config_path}"
    
    config = CGSimulationConfig.from_yaml(str(config_path))
    print(f"✓ 配置文件加载成功")
    print(f"  - 系统: {config.system_name}")
    print(f"  - 盒子: {config.box}")
    print(f"  - 拓扑: {config.topol}")
    print(f"  - 温度: {config.temperature} K")
    print(f"  - 离子强度: {config.ionic} M")
    
    # 验证 SLAB 拓扑
    assert config.topol == TopologyType.SLAB, f"期望 SLAB 拓扑，实际 {config.topol}"
    print(f"  ✓ SLAB 拓扑验证通过")
    
    return config

def test_config_validation(config):
    """测试 3.2: 配置验证"""
    print("\n============================================================")
    print("测试 3.2: SLAB 配置验证")
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
    assert comp.type == ComponentType.IDP, "应为 IDP 组件"
    print(f"  ✓ 组件验证通过")

def test_slab_params_auto_calculation():
    """测试 3.3: SLAB 参数自动计算"""
    print("\n============================================================")
    print("测试 3.3: SLAB 参数自动计算")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab.yaml"
    config = CGSimulationConfig.from_yaml(str(config_path))
    
    wrapper = CalvadosWrapper(config)
    
    # 生成配置并验证 SLAB 参数
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = os.path.join(tmpdir, config.system_name + "_CG")
        wrapper.write(output_dir, overwrite=True)
        
        # 读取生成的 config.yaml
        import yaml
        with open(os.path.join(output_dir, "config.yaml")) as f:
            calvados_config = yaml.safe_load(f)
        
        # 验证 SLAB 自动设置
        expected_slab_width = config.box[2] / 2  # z 方向的一半
        
        assert calvados_config.get('topol') == 'slab', "topol 应为 slab"
        assert calvados_config.get('slab_eq') == False, "slab_eq 应为 false（不使用 equilibration）"
        assert calvados_config.get('slab_width') == expected_slab_width, \
            f"slab_width 应为 {expected_slab_width}，实际 {calvados_config.get('slab_width')}"
        
        print(f"✓ SLAB 参数自动计算正确")
        print(f"  - topol: {calvados_config['topol']}")
        print(f"  - slab_eq: {calvados_config['slab_eq']}")
        print(f"  - slab_width: {calvados_config['slab_width']} (box[2]/2 = {expected_slab_width})")

def test_slab_simulation_run():
    """测试 3.4: 运行 SLAB 模拟"""
    print("\n============================================================")
    print("测试 3.4: 运行 SLAB 模拟")
    print("============================================================")
    
    config_path = Path(__file__).parent / "config_idp_slab.yaml"
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
    
    # 检查粒子数
    # 50 个分子 * 序列长度（约 100 残基）
    print(f"✓ SLAB 模拟完成!")
    print(f"  输出目录: {cg_output_dir}")
    print(f"  - final.pdb: ✓")
    print(f"  - trajectory.dcd: ✓")
    print(f"  - raw/: {len(os.listdir(raw_dir))} 个文件")

def main():
    """主测试函数"""
    print("=" * 60)
    print("SLAB 模拟测试 (config_idp_slab.yaml)")
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
        # 测试 3.1: 配置加载
        config = test_config_loading()
        
        # 测试 3.2: 配置验证
        test_config_validation(config)
        
        # 测试 3.3: SLAB 参数自动计算
        test_slab_params_auto_calculation()
        
        # 测试 3.4: 运行模拟
        test_slab_simulation_run()
        
        print("\n" + "=" * 60)
        print("🎉 所有 SLAB 测试通过!")
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
