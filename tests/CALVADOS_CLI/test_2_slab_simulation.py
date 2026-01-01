#!/usr/bin/env python3
"""
Test 2: CLI 单组分 SLAB 模拟测试

使用 CLI 命令测试单组分 IDP 在 SLAB 拓扑下的相分离模拟。
对比 test_3_slab_simulation.py 的脚本方式，这里使用 CLI 方式。
"""

import sys
import os
import subprocess
import shutil
from pathlib import Path


def test_cli_init_slab():
    """测试 2.1: 使用 CLI init 创建 SLAB 配置"""
    print("\n============================================================")
    print("测试 2.1: CLI init 创建 SLAB 配置")
    print("============================================================")

    # 使用临时目录
    test_dir = Path(__file__).parent / "test2_output"
    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 复制 FASTA 文件
        src_fasta = Path(__file__).parent.parent / "CALVADOS_IDP" / "TDP43_CTD.fasta"
        assert src_fasta.exists(), f"源 FASTA 文件不存在: {src_fasta}"
        shutil.copy2(src_fasta, test_dir / "TDP43_CTD.fasta")

        # 使用 CLI init 创建配置
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'init',
                'my_slab_simulation',
                '--type', 'idp',
                '--topol', 'slab',
                '--nmol', '50',
                '-o', str(test_dir)
            ],
            capture_output=True,
            text=True,
            timeout=30
        )

        assert result.returncode == 0, f"init 命令失败: {result.stderr}"
        assert 'Configuration template created' in result.stdout, "init 输出异常"
        print("✓ CLI init 成功")

        # 验证生成的配置文件
        config_file = test_dir / "my_slab_simulation.yaml"
        assert config_file.exists(), f"配置文件未生成: {config_file}"

        # 读取并验证配置
        from multiscale2.src import CGSimulationConfig, TopologyType
        config = CGSimulationConfig.from_yaml(str(config_file))

        assert config.system_name == 'my_slab_simulation', f"系统名错误: {config.system_name}"
        assert config.topol == TopologyType.SLAB, f"拓扑类型错误: {config.topol}"
        assert len(config.components) == 1, f"组件数错误: {len(config.components)}"
        assert config.components[0].name == 'protein_A', f"组件名错误: {config.components[0].name}"

        print(f"✓ 配置验证成功")
        print(f"  - 系统: {config.system_name}")
        print(f"  - 拓扑: {config.topol.value}")
        print(f"  - 组件: {config.components[0].name}")

        return True

    except subprocess.TimeoutExpired:
        print("✗ init 命令超时")
        return False
    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        os.chdir(original_cwd)


def test_cli_cg_simulation():
    """测试 2.2: 使用 CLI cg 运行 SLAB 模拟 (快速测试)"""
    print("\n============================================================")
    print("测试 2.2: CLI cg 运行 SLAB 模拟 (快速测试)")
    print("============================================================")

    # 使用 IDP 测试目录的配置文件
    src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
    test_dir = Path(__file__).parent / "test2_output"

    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 复制必要的输入文件
    src_config = src_dir / "config_idp_slab.yaml"
    src_fasta = src_dir / "TDP43_CTD.fasta"

    assert src_config.exists(), f"源配置文件不存在: {src_config}"
    assert src_fasta.exists(), f"源 FASTA 文件不存在: {src_fasta}"

    # 复制配置文件并修改路径
    import yaml
    with open(src_config, 'r') as f:
        config_data = yaml.safe_load(f)

    config_data['components'][0]['ffasta'] = 'TDP43_CTD.fasta'

    test_config = test_dir / "config_idp_slab.yaml"
    with open(test_config, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)

    shutil.copy2(src_fasta, test_dir / "TDP43_CTD.fasta")

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 运行 cg 命令
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'cg',
                '-f', str(test_config),
                '--overwrite'
            ],
            capture_output=True,
            text=True,
            timeout=300  # 5分钟超时
        )

        print(f"命令输出:\n{result.stdout}")
        if result.stderr:
            print(f"命令错误:\n{result.stderr}")

        assert result.returncode == 0, f"cg 命令失败: {result.stderr}"

        # 验证输出
        cg_output_dir = test_dir / f"{config_data['system_name']}_CG"
        assert cg_output_dir.exists(), f"输出目录不存在: {cg_output_dir}"

        # 检查关键文件
        assert (cg_output_dir / "final.pdb").exists(), "缺少 final.pdb"
        assert (cg_output_dir / "trajectory.xtc").exists(), "缺少 trajectory.xtc"
        assert (cg_output_dir / "simulation.log").exists(), "缺少 simulation.log"

        # 检查 raw 目录
        raw_dir = cg_output_dir / "raw"
        assert raw_dir.exists(), f"raw 目录不存在: {raw_dir}"

        print("✓ CLI cg 模拟成功")
        print(f"  输出目录: {cg_output_dir}")
        print(f"  - final.pdb: ✓")
        print(f"  - trajectory.xtc: ✓")
        print(f"  - raw/: {len(list(raw_dir.iterdir()))} 个文件")

        return True

    except subprocess.TimeoutExpired:
        print("✗ cg 命令超时")
        return False
    except AssertionError as e:
        print(f"✗ 验证失败: {e}")
        return False
    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        os.chdir(original_cwd)


def test_cli_cg_dry_run():
    """测试 2.3: CLI cg 干运行模式"""
    print("\n============================================================")
    print("测试 2.3: CLI cg 干运行模式")
    print("============================================================")

    src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
    test_dir = Path(__file__).parent / "test2_output"

    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 复制配置文件
    src_config = src_dir / "config_idp_slab.yaml"
    src_fasta = src_dir / "TDP43_CTD.fasta"

    shutil.copy2(src_config, test_dir / "config_idp_slab.yaml")
    shutil.copy2(src_fasta, test_dir / "TDP43_CTD.fasta")

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 运行干运行
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'cg',
                '-f', 'config_idp_slab.yaml',
                '--dry-run',
                '--overwrite'
            ],
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode != 0:
            print(f"命令输出:\n{result.stdout}")
            print(f"命令错误:\n{result.stderr}")
        assert result.returncode == 0, f"dry-run 命令失败"
        assert 'Dry run' in result.stdout or 'configuration files generated' in result.stdout, "dry-run 输出异常"

        print("✓ CLI cg 干运行成功")

        # 验证生成了配置文件（在 raw/ 目录下）
        output_dir = test_dir / "TDP43_CTD_slab_CG"
        raw_dir = output_dir / "raw"
        assert (raw_dir / "config.yaml").exists(), "dry-run 未生成 config.yaml"
        assert (raw_dir / "components.yaml").exists(), "dry-run 未生成 components.yaml"

        print(f"  生成了配置文件: ✓")
        print(f"  - config.yaml: ✓")
        print(f"  - components.yaml: ✓")

        return True

    except subprocess.TimeoutExpired:
        print("✗ dry-run 命令超时")
        return False
    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        os.chdir(original_cwd)


def test_cli_cg_with_custom_force_field():
    """测试 2.4: CLI cg 指定力场"""
    print("\n============================================================")
    print("测试 2.4: CLI cg 指定力场")
    print("============================================================")

    src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
    test_dir = Path(__file__).parent / "test2_output"

    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 复制配置文件
    src_config = src_dir / "config_idp_slab.yaml"
    src_fasta = src_dir / "TDP43_CTD.fasta"

    shutil.copy2(src_config, test_dir / "config_idp_slab.yaml")
    shutil.copy2(src_fasta, test_dir / "TDP43_CTD.fasta")

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 使用 calvados 力场运行（使用 --dry-run 避免长时间运行）
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'cg',
                '-f', 'config_idp_slab.yaml',
                '--force-field', 'calvados',
                '--dry-run',
                '--overwrite'
            ],
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode != 0:
            print(f"命令输出:\n{result.stdout}")
            print(f"命令错误:\n{result.stderr}")
        assert result.returncode == 0, f"cg --force-field 命令失败"
        assert 'calvados' in result.stdout.lower(), "力场输出异常"

        print("✓ CLI cg 力场指定成功")

        return True

    except subprocess.TimeoutExpired:
        print("✗ 命令超时")
        return False
    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        os.chdir(original_cwd)


def test_cli_cg_full_simulation():
    """测试 2.5: CLI cg 完整模拟 (使用配置中的步数)"""
    print("\n============================================================")
    print("测试 2.5: CLI cg 完整模拟")
    print("============================================================")

    src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
    test_dir = Path(__file__).parent / "test2_output_full"

    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 复制必要的输入文件
    src_config = src_dir / "config_idp_slab.yaml"
    src_fasta = src_dir / "TDP43_CTD.fasta"

    assert src_config.exists(), f"源配置文件不存在: {src_config}"
    assert src_fasta.exists(), f"源 FASTA 文件不存在: {src_fasta}"

    # 读取配置并修改步数
    import yaml
    with open(src_config, 'r') as f:
        config_data = yaml.safe_load(f)

    config_data['components'][0]['ffasta'] = 'TDP43_CTD.fasta'

    # 使用 config 中的步数进行测试
    config_steps = config_data.get('simulation', {}).get('steps') or config_data.get('steps', 20000)
    print(f"  配置步数: {config_steps}")

    test_config = test_dir / "config_idp_slab.yaml"
    with open(test_config, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)

    shutil.copy2(src_fasta, test_dir / "TDP43_CTD.fasta")

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 运行完整的 cg 命令（不使用 --dry-run）
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'cg',
                '-f', 'config_idp_slab.yaml',
                '--overwrite'
            ],
            capture_output=True,
            text=True,
            timeout=600  # 10分钟超时
        )

        print(f"命令输出:\n{result.stdout}")
        if result.stderr:
            print(f"命令错误:\n{result.stderr}")

        assert result.returncode == 0, f"cg 命令失败: {result.stderr}"

        # 验证输出
        cg_output_dir = test_dir / f"{config_data['system_name']}_CG"
        assert cg_output_dir.exists(), f"输出目录不存在: {cg_output_dir}"

        # 检查关键文件
        assert (cg_output_dir / "final.pdb").exists(), "缺少 final.pdb"
        assert (cg_output_dir / "trajectory.xtc").exists(), "缺少 trajectory.xtc"
        assert (cg_output_dir / "simulation.log").exists(), "缺少 simulation.log"

        # 检查 raw 目录
        raw_dir = cg_output_dir / "raw"
        assert raw_dir.exists(), f"raw 目录不存在: {raw_dir}"

        # 检查模拟日志是否存在（步数信息可能以不同格式存在）
        with open(cg_output_dir / "simulation.log", 'r') as f:
            log_content = f.read()
        assert len(log_content) > 0, "日志文件为空"

        print("✓ CLI cg 完整模拟成功")
        print(f"  输出目录: {cg_output_dir}")
        print(f"  - final.pdb: ✓")
        print(f"  - trajectory.xtc: ✓")
        print(f"  - simulation.log: ✓")
        print(f"  - raw/: {len(list(raw_dir.iterdir()))} 个文件")

        return True

    except subprocess.TimeoutExpired:
        print("✗ cg 命令超时")
        return False
    except AssertionError as e:
        print(f"✗ 验证失败: {e}")
        return False
    except Exception as e:
        print(f"✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        os.chdir(original_cwd)


def main():
    """主测试函数"""
    print("=" * 60)
    print("CALVADOS CLI 单组分 SLAB 模拟测试 (Test 2)")
    print("=" * 60)
    print("使用 CLI 命令替代脚本方式运行模拟")

    all_passed = True

    # 测试 2.1: init 创建配置
    try:
        if not test_cli_init_slab():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 2.1 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 2.2: cg 运行模拟
    try:
        if not test_cli_cg_simulation():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 2.2 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 2.3: dry-run 模式
    try:
        if not test_cli_cg_dry_run():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 2.3 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 2.4: 力场指定
    try:
        if not test_cli_cg_with_custom_force_field():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 2.4 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 2.5: 完整模拟
    try:
        if not test_cli_cg_full_simulation():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 2.5 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 所有 CLI 单组分 SLAB 测试通过!")
    else:
        print("⚠️ 部分测试失败")
    print("=" * 60)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
