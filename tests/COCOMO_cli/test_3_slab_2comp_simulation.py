#!/usr/bin/env python3
"""
Test 3: COCOMO CLI 双组分 SLAB 模拟测试（简化版）

仅测试配置加载和拓扑构建，跳过耗时的模拟步骤。
"""

import sys
import os
import shutil
from pathlib import Path


def test_cocomo_2comp_config():
    """测试 3.1: COCOMO 双组分配置"""
    print("\n============================================================")
    print("测试 3.1: COCOMO 双组分配置")
    print("============================================================")

    try:
        from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType, TopologyType
        from multiscale2.src.cg import CGSimulator

        # 加载配置文件
        config_file = Path(__file__).parent / "config_idp_slab_2comp_tdp43_tdp43.yaml"
        assert config_file.exists(), f"配置文件不存在: {config_file}"

        config = CGSimulationConfig.from_yaml(str(config_file))
        print(f"✓ 配置加载成功: {config.system_name}")

        # 验证配置
        errors = config.validate()
        if errors:
            print(f"⚠ 配置验证警告: {errors}")
        else:
            print(f"✓ 配置验证通过")

        print(f"  - 盒子: {config.box} nm")
        print(f"  - 温度: {config.temperature} K")
        print(f"  - 拓扑: {config.topol.value}")
        print(f"  - 组件数: {len(config.components)}")

        # 显示组件信息
        for comp in config.components:
            print(f"  - 组件: {comp.name}, 类型: {comp.type.value}, 分子数: {comp.nmol}")

        # 创建模拟器
        simulator = CGSimulator(config)
        print("✓ CGSimulator 创建成功")

        return True

    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_cocomo_2comp_topology():
    """测试 3.2: COCOMO 双组分拓扑信息"""
    print("\n============================================================")
    print("测试 3.2: COCOMO 双组分拓扑信息")
    print("============================================================")

    try:
        from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType, TopologyType
        from multiscale2.src.cg import CGSimulator

        # 加载配置
        config_file = Path(__file__).parent / "config_idp_slab_2comp_tdp43_tdp43.yaml"
        config = CGSimulationConfig.from_yaml(str(config_file))

        # 创建模拟器
        simulator = CGSimulator(config)

        # 获取全局序列
        global_seq = simulator.get_global_sequence()
        print(f"✓ 全局序列获取成功")
        print(f"  - 序列长度: {len(global_seq)} 残基")

        # 获取链 ID
        chain_ids = simulator.get_chain_ids()
        n_chains = len(set(chain_ids))
        print(f"  - 链数: {n_chains}")

        # 获取链标识符
        chain_identifiers = simulator.get_chain_identifiers()
        unique_chains = simulator.get_unique_chain_identifiers()
        print(f"  - 唯一链标识符数: {len(unique_chains)}")

        # 获取 folded domain 信息
        folded = simulator.get_folded_domains()
        n_folded = sum(folded)
        print(f"  - Folded domain 残基数: {n_folded}")

        # 验证链 ID 数量与序列长度匹配
        assert len(chain_ids) == len(global_seq), "链ID数量与序列长度不匹配"
        print("✓ 拓扑信息验证通过")

        return True

    except Exception as e:
        print(f"✗ 拓扑信息测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """主测试函数"""
    print("=" * 60)
    print("COCOMO CLI 双组分 SLAB 模拟测试 (Test 3) - 简化版")
    print("=" * 60)

    all_passed = True

    # 测试 3.1: 双组分配置
    try:
        if not test_cocomo_2comp_config():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 3.1 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 3.2: 拓扑信息
    try:
        if not test_cocomo_2comp_topology():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 3.2 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 注意：系统创建测试需要完整的模拟环境，跳过
    print("\n⚠️ 系统创建测试需要完整模拟环境，已跳过")

    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 所有 COCOMO 双组分测试通过!")
    else:
        print("⚠️ 部分测试失败")
    print("=" * 60)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
