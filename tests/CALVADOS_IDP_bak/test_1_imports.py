#!/usr/bin/env python3
"""
Test 1: Basic Import Tests
基础导入测试

测试模块是否能正确导入。
"""

import sys
from pathlib import Path


def test_import_cg_module():
    """测试 CG 模块导入"""
    print("=" * 60)
    print("测试 1.1: 导入 CG 模块")
    print("=" * 60)
    
    try:
        from multiscale2.src import (
            CGSimulationConfig,
            CGComponent,
            ComponentType,
            TopologyType,
            Platform,
            SimulationParams,
            SimulationResult,
            CGSimulator,
        )
        print("✓ CG 配置类导入成功")
        return True
    except Exception as e:
        print(f"✗ CG 模块导入失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_import_calvados_wrapper():
    """测试 CALVADOS wrapper 导入"""
    print("\n" + "=" * 60)
    print("测试 1.2: 导入 CALVADOS Wrapper")
    print("=" * 60)
    
    try:
        from multiscale2.src.calvados_wrapper import (
            CalvadosWrapper,
            run_calvados,
        )
        print("✓ CALVADOS wrapper 导入成功")
        return True
    except Exception as e:
        print(f"✗ CALVADOS wrapper 导入失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_import_ms2_calvados():
    """测试 ms2_calvados 导入"""
    print("\n" + "=" * 60)
    print("测试 1.3: 导入内化 ms2_calvados")
    print("=" * 60)
    
    try:
        from multiscale2.extern.ms2_calvados import Config, Components
        print("✓ ms2_calvados Config, Components 导入成功")
        
        from multiscale2.extern.ms2_calvados import sim
        print("✓ ms2_calvados sim 导入成功")
        return True
    except Exception as e:
        print(f"✗ ms2_calvados 导入失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_import_from_file():
    """测试从配置文件导入"""
    print("\n" + "=" * 60)
    print("测试 1.4: 从配置文件导入")
    print("=" * 60)
    
    try:
        # 配置文件路径
        config_path = Path(__file__).parent / "config_idp.yaml"
        
        if not config_path.exists():
            print(f"✗ 配置文件不存在: {config_path}")
            return False
        
        # 导入配置
        from multiscale2.src import CGSimulationConfig
        config = CGSimulationConfig.from_yaml(str(config_path))
        
        print(f"✓ 配置文件导入成功")
        print(f"  - 系统名称: {config.system_name}")
        print(f"  - 盒子尺寸: {config.box}")
        print(f"  - 温度: {config.temperature} K")
        print(f"  - 拓扑: {config.topol.value}")
        print(f"  - 组件数: {len(config.components)}")
        print(f"  - 总分子数: {config.total_molecules()}")
        
        # 验证配置
        errors = config.validate()
        if errors:
            print(f"✗ 配置验证失败: {errors}")
            return False
        
        print(f"  ✓ 配置验证通过")
        return True
        
    except Exception as e:
        print(f"✗ 配置文件导入失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """主测试函数"""
    print("\n" + "=" * 60)
    print("基础导入测试")
    print("=" * 60)
    print()
    
    results = []
    
    results.append(("CG 模块导入", test_import_cg_module()))
    results.append(("CALVADOS Wrapper 导入", test_import_calvados_wrapper()))
    results.append(("ms2_calvados 导入", test_import_ms2_calvados()))
    results.append(("配置文件导入", test_import_from_file()))
    
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
        print("🎉 所有基础测试通过!")
        return 0
    else:
        print("⚠️  部分测试失败")
        return 1


if __name__ == "__main__":
    sys.exit(main())

