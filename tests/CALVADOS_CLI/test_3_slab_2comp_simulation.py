#!/usr/bin/env python3
"""
Test 3: CLI 双组分 SLAB 模拟测试 (TDP43 + TDP43_CTD)

使用 CLI 命令测试 MDP + IDP 混合组件在 SLAB 拓扑下的相分离模拟。
对比 test_6_slab_2comp_tdp43_tdp43.py 的脚本方式，这里使用 CLI 方式。
"""

import sys
import os
import subprocess
import shutil
from pathlib import Path


def test_cli_init_mixed():
    """测试 3.1: 使用 CLI init 创建混合组件配置"""
    print("\n============================================================")
    print("测试 3.1: CLI init 创建混合组件配置")
    print("============================================================")

    test_dir = Path(__file__).parent / "test3_output"
    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 复制输入文件
        src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
        shutil.copy2(src_dir / "TDP43.pdb", test_dir / "protein_A.pdb")
        shutil.copy2(src_dir / "TDP43_CTD.fasta", test_dir / "protein_B.fasta")

        # 创建域定义文件
        domains_content = """
protein_A:
  - [3, 76]
  - [106, 176]
  - [192, 260]
  - [320, 334]
"""
        with open(test_dir / "protein_A_domains.yaml", 'w') as f:
            f.write(domains_content)

        # 使用 CLI init 创建混合配置
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'init',
                'mixed_simulation',
                '--type', 'mixed',
                '--topol', 'slab',
                '--nmol', '20',
                '-o', str(test_dir)
            ],
            capture_output=True,
            text=True,
            timeout=30
        )

        assert result.returncode == 0, f"init mixed 命令失败: {result.stderr}"
        assert 'Mixed IDP + MDP' in result.stdout, "init 输出异常"
        print("✓ CLI init 混合配置成功")

        # 验证配置文件
        config_file = test_dir / "mixed_simulation.yaml"
        assert config_file.exists(), f"配置文件未生成: {config_file}"

        from multiscale2.src import CGSimulationConfig, ComponentType, TopologyType
        config = CGSimulationConfig.from_yaml(str(config_file))

        assert config.topol == TopologyType.SLAB, f"拓扑类型错误"
        assert len(config.components) == 2, f"组件数错误: {len(config.components)}"

        # 验证组件类型
        comp_types = [c.type for c in config.components]
        assert ComponentType.IDP in comp_types, "缺少 IDP 组件"
        assert ComponentType.MDP in comp_types, "缺少 MDP 组件"

        print(f"✓ 混合配置验证成功")
        print(f"  - 拓扑: {config.topol.value}")
        print(f"  - 组件数: {len(config.components)}")

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


def test_cli_cg_2comp_simulation():
    """测试 3.2: CLI cg 运行双组分 SLAB 模拟"""
    print("\n============================================================")
    print("测试 3.2: CLI cg 运行双组分 SLAB 模拟 (TDP43 + TDP43_CTD)")
    print("============================================================")

    src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
    test_dir = Path(__file__).parent / "test3_output"

    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 复制输入文件
    shutil.copy2(src_dir / "TDP43.pdb", test_dir / "TDP43.pdb")
    shutil.copy2(src_dir / "TDP43_CTD.fasta", test_dir / "TDP43_CTD.fasta")

    # 读取并修改配置
    import yaml
    with open(src_dir / "config_idp_slab_2comp_tdp43_tdp43.yaml", 'r') as f:
        config_data = yaml.safe_load(f)

    # 修正路径
    config_data['components'][0]['fpdb'] = 'TDP43.pdb'
    config_data['components'][0]['ffasta'] = 'TDP43_CTD.fasta'
    config_data['components'][1]['ffasta'] = 'TDP43_CTD.fasta'

    test_config = test_dir / "config_2comp.yaml"
    with open(test_config, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)

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
            timeout=600  # 10分钟超时，双组分更复杂
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

        # 验证 simulation.log 包含双组分信息
        with open(cg_output_dir / "simulation.log", 'r') as f:
            log_content = f.read()
        assert 'TDP43' in log_content, "日志中缺少 TDP43 组件"
        assert 'TDP43_CTD' in log_content, "日志中缺少 TDP43_CTD 组件"

        print("✓ CLI cg 双组分模拟成功")
        print(f"  输出目录: {cg_output_dir}")
        print(f"  - final.pdb: ✓")
        print(f"  - trajectory.xtc: ✓")
        print(f"  - simulation.log 包含双组分: ✓")
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


def test_cli_cg_2comp_dry_run():
    """测试 3.3: CLI cg 双组分干运行验证配置生成"""
    print("\n============================================================")
    print("测试 3.3: CLI cg 双组分干运行验证配置生成")
    print("============================================================")

    src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
    test_dir = Path(__file__).parent / "test3_output"

    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 复制输入文件
    shutil.copy2(src_dir / "TDP43.pdb", test_dir / "TDP43.pdb")
    shutil.copy2(src_dir / "TDP43_CTD.fasta", test_dir / "TDP43_CTD.fasta")

    # 读取并修改配置
    import yaml
    with open(src_dir / "config_idp_slab_2comp_tdp43_tdp43.yaml", 'r') as f:
        config_data = yaml.safe_load(f)

    config_data['components'][0]['fpdb'] = 'TDP43.pdb'
    config_data['components'][0]['ffasta'] = 'TDP43_CTD.fasta'
    config_data['components'][1]['ffasta'] = 'TDP43_CTD.fasta'

    test_config = test_dir / "config_2comp.yaml"
    with open(test_config, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 运行干运行
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'cg',
                '-f', str(test_config),
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
        print("✓ CLI cg 双组分干运行成功")

        # 验证生成了配置文件（在 raw/ 目录下）
        output_dir = test_dir / f"{config_data['system_name']}_CG"
        raw_dir = output_dir / "raw"
        assert (raw_dir / "config.yaml").exists(), "dry-run 未生成 config.yaml"
        assert (raw_dir / "components.yaml").exists(), "dry-run 未生成 components.yaml"

        # 验证 components.yaml 包含两个组件
        with open(raw_dir / "components.yaml", 'r') as f:
            components_content = f.read()
        assert 'TDP43' in components_content, "components.yaml 缺少 TDP43"
        assert 'TDP43_CTD' in components_content, "components.yaml 缺少 TDP43_CTD"

        print(f"  双组分配置文件验证: ✓")

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


def test_cli_inline_fdomains():
    """测试 3.4: CLI cg 内联 fdomains 处理"""
    print("\n============================================================")
    print("测试 3.4: CLI cg 内联 fdomains 处理")
    print("============================================================")

    test_dir = Path(__file__).parent / "test3_output"

    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 复制输入文件
    src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
    shutil.copy2(src_dir / "TDP43.pdb", test_dir / "TDP43.pdb")
    shutil.copy2(src_dir / "TDP43_CTD.fasta", test_dir / "TDP43_CTD.fasta")

    # 创建包含内联 fdomains 的配置
    config_content = """
system_name: TDP43_CTD_TDP43_SLAB
box: [40, 40, 80]
temperature: 293.0
ionic: 0.15
steps: 20000
report_interval: 2000
timestep: 0.01
topol: slab

components:
  - name: TDP43
    type: mdp
    nmol: 50
    fpdb: TDP43.pdb
    ffasta: TDP43_CTD.fasta
    restraint: true
    fdomains: |
      TDP43:
        - [3, 76]
        - [106, 176]
        - [192, 260]
        - [320, 334]
        
  - name: TDP43_CTD
    type: idp
    nmol: 50
    ffasta: TDP43_CTD.fasta
    charge_termini: both
    restraint: false
"""
    test_config = test_dir / "config_inline.yaml"
    with open(test_config, 'w') as f:
        f.write(config_content)

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 运行干运行验证内联 fdomains 被处理
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'cg',
                '-f', str(test_config),
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

        # 验证生成了域文件（在 raw/ 目录下）
        output_dir = test_dir / "TDP43_CTD_TDP43_SLAB_CG"
        raw_dir = output_dir / "raw"
        tdp43_domains = raw_dir / "TDP43_domains.yaml"
        assert tdp43_domains.exists(), f"内联 fdomains 未被处理为文件: {tdp43_domains}"

        # 验证域文件内容
        with open(tdp43_domains, 'r') as f:
            domains_content = f.read()
        assert 'TDP43' in domains_content, "域文件内容异常"
        assert '3' in domains_content and '76' in domains_content, "域定义缺失"

        print("✓ CLI cg 内联 fdomains 处理成功")
        print(f"  生成文件: {tdp43_domains.name}")

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


def test_cli_cg_2comp_full_simulation():
    """测试 3.5: CLI cg 双组分完整模拟 (使用配置中的步数)"""
    print("\n============================================================")
    print("测试 3.5: CLI cg 双组分完整模拟")
    print("============================================================")

    src_dir = Path(__file__).parent.parent / "CALVADOS_IDP"
    test_dir = Path(__file__).parent / "test3_output_full"

    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir(parents=True, exist_ok=True)

    # 复制输入文件
    shutil.copy2(src_dir / "TDP43.pdb", test_dir / "TDP43.pdb")
    shutil.copy2(src_dir / "TDP43_CTD.fasta", test_dir / "TDP43_CTD.fasta")

    # 读取配置
    import yaml
    with open(src_dir / "config_idp_slab_2comp_tdp43_tdp43.yaml", 'r') as f:
        config_data = yaml.safe_load(f)

    # 修正路径
    config_data['components'][0]['fpdb'] = 'TDP43.pdb'
    config_data['components'][0]['ffasta'] = 'TDP43_CTD.fasta'
    config_data['components'][1]['ffasta'] = 'TDP43_CTD.fasta'

    # 获取配置中的步数
    config_steps = config_data.get('simulation', {}).get('steps') or config_data.get('steps', 20000)
    print(f"  配置步数: {config_steps}")

    test_config = test_dir / "config_2comp.yaml"
    with open(test_config, 'w') as f:
        yaml.dump(config_data, f, default_flow_style=False)

    # 切换到测试目录
    original_cwd = os.getcwd()
    os.chdir(test_dir)

    try:
        # 运行完整的 cg 命令
        result = subprocess.run(
            [
                sys.executable, '-m', 'multiscale2.cli', 'cg',
                '-f', 'config_2comp.yaml',
                '--overwrite'
            ],
            capture_output=True,
            text=True,
            timeout=900  # 15分钟超时，双组分更复杂
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

        # 验证 simulation.log 包含双组分信息
        with open(cg_output_dir / "simulation.log", 'r') as f:
            log_content = f.read()
        assert 'TDP43' in log_content, "日志中缺少 TDP43 组件"
        assert 'TDP43_CTD' in log_content, "日志中缺少 TDP43_CTD 组件"
        assert len(log_content) > 0, "日志文件为空"

        print("✓ CLI cg 双组分完整模拟成功")
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
    print("CALVADOS CLI 双组分 SLAB 模拟测试 (Test 3)")
    print("=" * 60)
    print("使用 CLI 命令替代脚本方式运行 MDP + IDP 混合模拟")

    all_passed = True

    # 测试 3.1: init 创建混合配置
    try:
        if not test_cli_init_mixed():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 3.1 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 3.2: cg 运行双组分模拟
    try:
        if not test_cli_cg_2comp_simulation():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 3.2 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 3.3: 双组分干运行
    try:
        if not test_cli_cg_2comp_dry_run():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 3.3 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 3.4: 内联 fdomains 处理
    try:
        if not test_cli_inline_fdomains():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 3.4 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    # 测试 3.5: 双组分完整模拟
    try:
        if not test_cli_cg_2comp_full_simulation():
            all_passed = False
    except Exception as e:
        print(f"✗ 测试 3.5 异常: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False

    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 所有 CLI 双组分 SLAB 测试通过!")
    else:
        print("⚠️ 部分测试失败")
    print("=" * 60)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

