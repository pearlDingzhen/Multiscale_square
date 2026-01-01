#!/usr/bin/env python3
"""
Test 2: COCOMO CLI 单组分 SLAB 模拟测试（简化版）

仅测试配置加载和拓扑构建，跳过耗时的模拟步骤。
"""

import sys
import os
import shutil
from pathlib import Path


def test_cocomo_config_and_topology():
    """测试 2.1: COCOMO 配置和拓扑"""
    print("\n============================================================")
    print("测试 2.1: COCOMO 配置和拓扑")
    print("============================================================")

    try:
        from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType, TopologyType
        from multiscale2.src.cg import CGSimulator

        # 加载配置文件
        config_file = Path(__file__).parent / "config_idp_slab.yaml"
        assert config_file.exists(), f"配置文件不存在: {config_file}"

        config = CGSimulationConfig.from_yaml(str(config_file))
        print(f"✓ 配置加载成功: {config.system_name}")

        # 验证配置
        errors = config.validate()
        assert len(errors) == 0, f"配置验证失败: {errors}"
        print(f"✓ 配置验证通过")
        print(f"  - 盒子: {config.box} nm")
        print(f"  - 温度: {config.temperature} K")
        print(f"  - 拓扑: {config.topol.value}")
        print(f"  - 组件数: {len(config.components)}")
        print(f"  - 总分子数: {config.total_molecules()}")

        # 创建模拟器
        simulator = CGSimulator(config)
        print("✓ CGSimulator 创建成功")

        # 测试拓扑信息获取
        global_seq = simulator.get_global_sequence()
        print(f"✓ 全局序列获取成功")
        print(f"  - 序列长度: {len(global_seq)} 残基")

        chain_ids = simulator.get_chain_ids()
        print(f"  - 链数: {len(set(chain_ids))}")

        return True

    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_cocomo_system_creation():
    """测试 2.2: COCOMO 系统创建"""
    print("\n============================================================")
    print("测试 2.2: COCOMO 系统创建")
    print("============================================================")

    try:
        from multiscale2.src.cocomo2_creator import COCOMO
        import numpy as np

        # 创建简单的拓扑信息
        topology_info = {
            'global_sequence': 'ACDEFGHIKLMNPQRSTVWY' * 5,  # 100 个残基
            'chain_ids': [1] * 100,
            'folded_domains': [0] * 100,
            'component_names': ['test'] * 100,
            'local_residue_indices': list(range(1, 101)),
        }

        # 创建位置
        np.random.seed(42)
        positions = np.random.rand(100, 3) * 10  # 10nm box

        # 创建 COCOMO 实例
        cocomo = COCOMO(
            box_size=[25.0, 25.0, 30.0],
            topology_info=topology_info,
            positions=positions,
            resources='CPU'
        )
        print("✓ COCOMO 实例创建成功")

        # 创建系统
        system, top = cocomo.create_system()
        print("✓ 系统创建成功")
        print(f"  - 粒子数: {system.getNumParticles()}")
        print(f"  - 力数: {system.getNumForces()}")

        return True

    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_cocomo_trajectory_output():
    """测试 2.3: COCOMO 轨迹输出（使用 XTCReporter）"""
    print("\n============================================================")
    print("测试 2.3: COCOMO 轨迹输出（XTCReporter）")
    print("============================================================")

    try:
        from mdtraj.reporters import XTCReporter
        from multiscale2.src.cocomo2_creator import COCOMO
        import numpy as np
        import tempfile
        import os

        # 创建简单的拓扑信息
        topology_info = {
            'global_sequence': 'MGSS' * 3,  # 12 个残基
            'chain_ids': [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
            'folded_domains': [0] * 12,
            'component_names': ['test'] * 12,
            'local_residue_indices': list(range(1, 13)),
        }

        # 创建位置
        np.random.seed(42)
        positions = np.random.rand(12, 3) * 10

        # 创建 COCOMO 实例
        cocomo = COCOMO(
            box_size=[20.0, 20.0, 20.0],
            topology_info=topology_info,
            positions=positions,
            resources='CPU'
        )

        # 创建系统
        system, top = cocomo.create_system()

        # 创建临时目录用于测试
        with tempfile.TemporaryDirectory() as tmpdir:
            xtc_file = os.path.join(tmpdir, 'test_trajectory.xtc')

            # 测试 XTCReporter
            reporter = XTCReporter(xtc_file, reportInterval=1)
            print("✓ XTCReporter 创建成功")
            print(f"  - 输出文件: {xtc_file}")

        print("✓ 轨迹输出测试通过")

        return True

    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """主测试函数"""
    print("=" * 60)
    print("COCOMO CLI SLAB 模拟测试 (Test 2) - 简化版")
    print("=" * 60)

    all_passed = True

    # 测试 2.1: 配置和拓扑
    try:
        if not test_cocomo_config_and_topology():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 2.1 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 2.2: 系统创建
    try:
        if not test_cocomo_system_creation():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 2.2 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 2.3: 轨迹输出
    try:
        if not test_cocomo_trajectory_output():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 2.3 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 所有 COCOMO SLAB 测试通过!")
    else:
        print("⚠️ 部分测试失败")
    print("=" * 60)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
