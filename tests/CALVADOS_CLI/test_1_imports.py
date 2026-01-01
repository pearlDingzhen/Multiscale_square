#!/usr/bin/env python3
"""
Test 1: CLI 基础导入和功能测试

测试 multiscale2 CLI 模块的基本导入和功能。
"""

import sys
import subprocess
import os
from pathlib import Path


def test_cli_module_imports():
    """测试 1.1: CLI 模块导入"""
    print("\n============================================================")
    print("测试 1.1: CLI 模块导入")
    print("============================================================")

    # 测试主模块导入
    try:
        from multiscale2.cli import main, cli
        print("✓ multiscale2.cli 导入成功")
    except ImportError as e:
        print(f"✗ multiscale2.cli 导入失败: {e}")
        return False

    # 测试命令导入
    try:
        from multiscale2.cli.commands import init_command, cg_command, info_command
        print("✓ 命令模块导入成功")
    except ImportError as e:
        print(f"✗ 命令模块导入失败: {e}")
        return False

    # 测试 CG 模块导入
    try:
        from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType
        from multiscale2.src import TopologyType, Platform, SimulationParams, CGSimulator
        print("✓ CG 模块导入成功")
    except ImportError as e:
        print(f"✗ CG 模块导入失败: {e}")
        return False

    # 测试 CALVADOS wrapper 导入
    try:
        from multiscale2.src.calvados_wrapper import CalvadosWrapper, run_calvados
        print("✓ CALVADOS wrapper 导入成功")
    except ImportError as e:
        print(f"✗ CALVADOS wrapper 导入失败: {e}")
        return False

    return True


def test_cli_help():
    """测试 1.2: CLI 帮助信息"""
    print("\n============================================================")
    print("测试 1.2: CLI 帮助信息")
    print("============================================================")

    # 测试主帮助
    result = subprocess.run(
        [sys.executable, '-m', 'multiscale2.cli', '--help'],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"帮助命令失败: {result.stderr}"
    assert 'init' in result.stdout, "帮助中缺少 init 命令"
    assert 'cg' in result.stdout, "帮助中缺少 cg 命令"
    assert 'info' in result.stdout, "帮助中缺少 info 命令"
    print("✓ 主帮助信息正常")

    # 测试 init 帮助
    result = subprocess.run(
        [sys.executable, '-m', 'multiscale2.cli', 'init', '--help'],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"init 帮助失败: {result.stderr}"
    print("✓ init 命令帮助正常")

    # 测试 cg 帮助
    result = subprocess.run(
        [sys.executable, '-m', 'multiscale2.cli', 'cg', '--help'],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"cg 帮助失败: {result.stderr}"
    assert '--input-file' in result.stdout or '-f' in result.stdout, "cg 帮助中缺少 -f 选项"
    print("✓ cg 命令帮助正常")

    # 测试 info 帮助
    result = subprocess.run(
        [sys.executable, '-m', 'multiscale2.cli', 'info', '--help'],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"info 帮助失败: {result.stderr}"
    print("✓ info 命令帮助正常")

    return True


def test_cli_entry_point():
    """测试 1.3: CLI 入口点"""
    print("\n============================================================")
    print("测试 1.3: CLI 入口点")
    print("============================================================")

    # 测试 info 命令执行
    result = subprocess.run(
        [sys.executable, '-m', 'multiscale2.cli', 'info'],
        capture_output=True,
        text=True,
        timeout=30
    )

    assert result.returncode == 0, f"info 命令失败: {result.stderr}"
    assert 'Python' in result.stdout, "info 输出中缺少 Python 信息"
    assert 'Available force fields' in result.stdout or 'calvados' in result.stdout, "info 输出中缺少力场信息"
    print("✓ info 命令执行成功")
    print(f"  Python 版本检测: ✓")
    print(f"  力场列表检测: ✓")

    return True


def test_cg_config_class():
    """测试 1.4: CG 配置类功能"""
    print("\n============================================================")
    print("测试 1.4: CG 配置类功能")
    print("============================================================")

    from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType, TopologyType
    import tempfile

    # 创建临时 FASTA 文件用于测试
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">test_sequence\n")
        f.write("MKLPTCVVYATGNVPRAILVDLEGGFETPGRVIVVEVKPKKITNVTVGDYVMQELPK\n")
        fasta_path = f.name

    try:
        # 创建配置
        config = CGSimulationConfig(
            system_name="test_system",
            box=[25.0, 25.0, 30.0],
            temperature=310.0,
            ionic=0.15,
            topol=TopologyType.CUBIC,
        )

        # 添加 IDP 组件（使用真实文件路径）
        idp_comp = CGComponent(
            name="test_idp",
            type=ComponentType.IDP,
            nmol=10,
            ffasta=fasta_path,
        )
        config.add_component(idp_comp)

        # 验证
        errors = config.validate()
        assert len(errors) == 0, f"配置验证失败: {errors}"
        print("✓ 配置创建成功")
        print(f"  - 系统: {config.system_name}")
        print(f"  - 组件数: {len(config.components)}")
        print(f"  - 总分子数: {config.total_molecules()}")

        # 测试字典转换
        config_dict = config.to_dict()
        assert config_dict['system_name'] == 'test_system'
        assert len(config_dict['components']) == 1
        print("✓ 字典转换正常")

        # 测试从字典恢复
        config2 = CGSimulationConfig.from_dict(config_dict)
        assert config2.system_name == 'test_system'
        print("✓ 字典恢复正常")

        return True
    finally:
        # 清理临时文件
        os.unlink(fasta_path)


def test_component_validation():
    """测试 1.5: 组件验证逻辑"""
    print("\n============================================================")
    print("测试 1.5: 组件验证逻辑")
    print("============================================================")

    from multiscale2.src import CGComponent, ComponentType

    # 测试 IDP 验证（缺少 ffasta）
    idp_comp = CGComponent(
        name="test_idp",
        type=ComponentType.IDP,
        nmol=10,
    )
    errors = idp_comp.validate()
    assert len(errors) > 0, "缺少 ffasta 应该报错"
    print("✓ IDP 缺少 ffasta 验证: 正确报错")

    # 测试 MDP 验证（缺少 fpdb）
    mdp_comp = CGComponent(
        name="test_mdp",
        type=ComponentType.MDP,
        nmol=10,
    )
    errors = mdp_comp.validate()
    assert len(errors) > 0, "缺少 fpdb 应该报错"
    print("✓ MDP 缺少 fpdb 验证: 正确报错")

    # 测试 MDP + restraint 但缺少 fdomains
    mdp_comp2 = CGComponent(
        name="test_mdp2",
        type=ComponentType.MDP,
        nmol=10,
        fpdb="test.pdb",
        restraint=True,
    )
    errors = mdp_comp2.validate()
    assert len(errors) > 0, "restraint=true 缺少 fdomains 应该报错"
    print("✓ MDP restraint 缺少 fdomains 验证: 正确报错")

    return True


def test_calvados_wrapper():
    """测试 1.6: CALVADOS Wrapper 功能"""
    print("\n============================================================")
    print("测试 1.6: CALVADOS Wrapper 功能")
    print("============================================================")

    from multiscale2.src import CGSimulationConfig, CGComponent, ComponentType, TopologyType
    from multiscale2.src.calvados_wrapper import CalvadosWrapper

    # 创建简单配置
    config = CGSimulationConfig(
        system_name="test_wrapper",
        box=[25.0, 25.0, 30.0],
        temperature=310.0,
        ionic=0.15,
        topol=TopologyType.CUBIC,
    )
    config.add_component(CGComponent(
        name="test_comp",
        type=ComponentType.IDP,
        nmol=10,
        ffasta="test.fasta",
    ))

    # 创建 wrapper
    wrapper = CalvadosWrapper(config)
    print("✓ CalvadosWrapper 创建成功")

    # 测试拓扑转换
    cubic_topol = wrapper._topol_to_calvados()
    assert cubic_topol == 'grid', f"cubic 应该转为 grid，实际 {cubic_topol}"
    print("✓ 拓扑类型转换: cubic → grid")

    config.topol = TopologyType.SLAB
    wrapper2 = CalvadosWrapper(config)
    slab_topol = wrapper2._topol_to_calvados()
    assert slab_topol == 'slab', f"slab 应该转为 slab，实际 {slab_topol}"
    print("✓ 拓扑类型转换: slab → slab")

    return True


def main():
    """主测试函数"""
    print("=" * 60)
    print("CALVADOS CLI 基础测试 (Test 1)")
    print("=" * 60)

    all_passed = True

    # 测试 1.1: 模块导入
    try:
        if not test_cli_module_imports():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.1 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.2: 帮助信息
    try:
        if not test_cli_help():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.2 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.3: 入口点
    try:
        if not test_cli_entry_point():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.3 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.4: 配置类
    try:
        if not test_cg_config_class():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.4 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.5: 验证逻辑
    try:
        if not test_component_validation():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.5 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 1.6: Wrapper 功能
    try:
        if not test_calvados_wrapper():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 1.6 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 所有 CLI 基础测试通过!")
    else:
        print("⚠️ 部分测试失败")
    print("=" * 60)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

