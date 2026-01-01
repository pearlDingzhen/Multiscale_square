#!/usr/bin/env python3
"""
Test 1: COCOMO CLI 基础导入和功能测试

测试 COCOMO CLI 模块的基本导入和功能。
"""

import sys
import subprocess
import os
from pathlib import Path


def test_cocomo_module_imports():
    """测试 1.1: COCOMO 模块导入"""
    print("\n============================================================")
    print("测试 1.1: COCOMO 模块导入")
    print("============================================================")

    # 测试 COCOMO 类导入
    try:
        from multiscale2.src.cocomo2_creator import COCOMO
        print("✓ COCOMO 类导入成功")
    except ImportError as e:
        print(f"✗ COCOMO 类导入失败: {e}")
        return False

    # 测试 CG 模块导入
    try:
        from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType
        from multiscale2.src import TopologyType, SimulationParams, CGSimulator
        print("✓ CG 模块导入成功")
    except ImportError as e:
        print(f"✗ CG 模块导入失败: {e}")
        return False

    # 测试 OpenMM 相关导入
    try:
        from openmm import Platform, LangevinIntegrator, XmlSerializer
        from openmm.app import Simulation, StateDataReporter
        print("✓ OpenMM 模块导入成功")
    except ImportError as e:
        print(f"✗ OpenMM 模块导入失败: {e}")
        return False

    # 测试 mdtraj reporters
    try:
        from mdtraj.reporters import XTCReporter
        print("✓ mdtraj XTCReporter 导入成功")
    except ImportError as e:
        print(f"✗ mdtraj XTCReporter 导入失败: {e}")
        return False

    return True


def test_cocomo_class_instantiation():
    """测试 1.2: COCOMO 类实例化"""
    print("\n============================================================")
    print("测试 1.2: COCOMO 类实例化")
    print("============================================================")

    from multiscale2.src.cocomo2_creator import COCOMO
    import numpy as np

    # 测试基本实例化
    try:
        # 创建简单的拓扑信息
        topology_info = {
            'global_sequence': 'ACDEFGHIKLMNPQRSTVWY' * 10,  # 200 个残基
            'chain_ids': [1] * 200,
            'folded_domains': [0] * 200,
            'component_names': ['test'] * 200,
            'local_residue_indices': list(range(1, 201)),
        }

        # 创建位置 (简单测试)
        positions = np.random.rand(200, 3) * 10  # 10nm box

        # 创建 COCOMO 实例
        cocomo = COCOMO(
            box_size=[25.0, 25.0, 30.0],
            topology_info=topology_info,
            positions=positions,
            resources='CPU'  # 测试用 CPU
        )

        print("✓ COCOMO 实例创建成功")
        print(f"  - 盒子大小: {cocomo.box_size}")
        print(f"  - 残基数: {len(topology_info['global_sequence'])}")

        return True
    except Exception as e:
        print(f"✗ COCOMO 实例化失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_cocomo_topology_conversion():
    """测试 1.3: COCOMO 拓扑转换"""
    print("\n============================================================")
    print("测试 1.3: COCOMO 拓扑转换")
    print("============================================================")

    from multiscale2.src.cocomo2_creator import COCOMO
    import numpy as np

    try:
        topology_info = {
            'global_sequence': 'MSEQ' * 5,
            'chain_ids': [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
            'folded_domains': [0] * 10,
            'component_names': ['test'] * 10,
            'local_residue_indices': list(range(1, 11)),
        }
        positions = np.random.rand(10, 3) * 10

        cocomo = COCOMO(
            box_size=[25.0, 25.0, 30.0],
            topology_info=topology_info,
            positions=positions,
            resources='CPU'
        )

        # 测试拓扑构建
        top = cocomo._build_topology()
        print("✓ 拓扑构建成功")
        print(f"  - 残基数: {top.getNumAtoms()}")
        print(f"  - 链数: {top.getNumChains()}")

        # 测试盒子向量创建
        box_vectors = cocomo._create_box_vectors()
        print("✓ 盒子向量创建成功")
        print(f"  - 向量数: {len(box_vectors)}")

        return True
    except Exception as e:
        print(f"✗ 拓扑转换测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_cg_simulator_cocomo():
    """测试 1.4: CGSimulator COCOMO 功能"""
    print("\n============================================================")
    print("测试 1.4: CGSimulator COCOMO 功能")
    print("============================================================")

    from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType, TopologyType
    from multiscale2.src.cg import CGSimulator
    import tempfile
    import os

    # 创建临时 FASTA 文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">test_seq\n")
        f.write("MKLPTCVVYATGNVPRAILVDLEGGFETPGRVIVVEVKPKKITNVTVGDYVMQELPK\n")
        fasta_path = f.name

    try:
        # 创建配置
        config = CGSimulationConfig(
            system_name="test_cocomo",
            box=[25.0, 25.0, 30.0],
            temperature=310.0,
            ionic=0.15,
            topol=TopologyType.SLAB,
        )
        config.add_component(CGComponent(
            name="test_idp",
            type=ComponentType.IDP,
            nmol=5,  # 少量分子用于测试
            ffasta=fasta_path,
        ))

        # 创建模拟器
        simulator = CGSimulator(config)
        print("✓ CGSimulator 创建成功")
        print(f"  - 系统: {config.system_name}")
        print(f"  - 总分子数: {config.total_molecules()}")

        # 测试获取组成信息
        composition = simulator.get_composition()
        assert len(composition) == 1
        print(f"  - 组件: {composition[0]['name']}")
        print(f"  - 序列长度: {composition[0]['nres']}")

        return True
    except Exception as e:
        print(f"✗ CGSimulator COCOMO 功能测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        os.unlink(fasta_path)


def test_cli_help_cocomo():
    """测试 1.5: CLI 帮助信息"""
    print("\n============================================================")
    print("测试 1.5: CLI 帮助信息")
    print("============================================================")

    # 测试主帮助
    result = subprocess.run(
        [sys.executable, '-m', 'multiscale2.cli', '--help'],
        capture_output=True,
        text=True,
        timeout=30
    )

    assert result.returncode == 0, f"帮助命令失败: {result.stderr}"
    print("✓ 主帮助信息正常")

    # 测试 info 命令包含 COCOMO
    result = subprocess.run(
        [sys.executable, '-m', 'multiscale2.cli', 'info'],
        capture_output=True,
        text=True,
        timeout=30
    )

    assert result.returncode == 0, f"info 命令失败: {result.stderr}"
    if 'cocomo' in result.stdout.lower() or 'COCOMO' in result.stdout:
        print("✓ info 命令包含 COCOMO 信息")
    else:
        print("⚠ info 命令可能未包含 COCOMO 信息")

    return True


def main():
    """主测试函数"""
    print("=" * 60)
    print("COCOMO CLI 基础测试 (Test 1)")
    print("=" * 60)

    all_passed = True

    # 测试 1.1: 模块导入
    try:
        if not test_cocomo_module_imports():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.1 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.2: COCOMO 实例化
    try:
        if not test_cocomo_class_instantiation():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.2 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.3: 拓扑转换
    try:
        if not test_cocomo_topology_conversion():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.3 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.4: CGSimulator COCOMO
    try:
        if not test_cg_simulator_cocomo():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.4 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.5: CLI 帮助
    try:
        if not test_cli_help_cocomo():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.5 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 所有 COCOMO CLI 基础测试通过!")
    else:
        print("⚠️ 部分测试失败")
    print("=" * 60)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

