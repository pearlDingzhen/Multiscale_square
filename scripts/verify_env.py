#!/usr/bin/env python3
"""
Multiscale² 环境验证脚本
"""

import sys


def check_versions():
    """检查所有核心包的版本"""
    print("=" * 60)
    print("Multiscale² 环境验证")
    print("=" * 60)

    # 核心包列表
    packages = [
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("scipy", "scipy"),
        ("openmm", "openmm"),
        ("mdtraj", "mdtraj"),
        ("mdanalysis", "MDAnalysis"),
        ("biopython", "biopython"),
        ("numba", "numba"),
        ("torch", "torch"),
        ("dgl", "dgl"),
        ("e3nn", "e3nn"),
        ("se3_transformer", "se3_transformer"),
        ("networkx", "networkx"),
        ("matplotlib", "matplotlib"),
    ]

    print("\n[核心包版本]")
    print("-" * 40)
    all_ok = True
    for name, import_name in packages:
        try:
            mod = __import__(import_name)
            version = getattr(mod, "__version__", "unknown")
            print(f"  ✓ {name:18} : {version}")
        except ImportError as e:
            print(f"  ✗ {name:18} : 未安装")
            all_ok = False

    # GPU 检查
    print("\n[GPU 状态]")
    print("-" * 40)
    try:
        import torch

        print(f"  PyTorch 版本    : {torch.__version__}")
        print(f"  CUDA 可用       : {torch.cuda.is_available()}")

        if torch.cuda.is_available():
            print(f"  CUDA 版本       : {torch.version.cuda}")
            print(f"  GPU 数量        : {torch.cuda.device_count()}")

            for i in range(torch.cuda.device_count()):
                gpu_name = torch.cuda.get_device_name(i)
                print(f"    GPU {i}         : {gpu_name}")
        else:
            print("  ⚠️  CUDA 不可用，将使用 CPU")
    except ImportError:
        print("  ✗ PyTorch 未安装")
        all_ok = False

    # Python 环境
    print("\n[Python 环境]")
    print("-" * 40)
    print(f"  Python 版本     : {sys.version}")
    print(f"  执行路径        : {sys.executable}")

    # 总结
    print("\n" + "=" * 60)
    if all_ok:
        print("✅ 环境检查通过!")
    else:
        print("⚠️  部分包未安装，请检查依赖")
    print("=" * 60)

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(check_versions())





