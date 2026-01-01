#!/bin/bash
# ===========================================
# Multiscale² 安装脚本
# ===========================================

set -e  # 遇到错误立即退出

echo "==========================================="
echo "Multiscale² 环境安装脚本"
echo "==========================================="

# 检查 Python 版本
PYTHON_VERSION=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
echo "[1/6] 当前 Python 版本: ${PYTHON_VERSION}"

if [ "$PYTHON_VERSION" != "3.9" ]; then
    echo "⚠️  警告: 推荐使用 Python 3.9，当前版本为 ${PYTHON_VERSION}"
    echo "继续安装可能遇到兼容性问题。"
    read -p "是否继续? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "安装已取消。"
        exit 1
    fi
fi

# 创建 conda 环境 (如果不存在)
if [ -z "$CONDA_DEFAULT_ENV" ] || [ "$CONDA_DEFAULT_ENV" != "multiscale2" ]; then
    echo "[2/6] 创建 conda 环境 multiscale2..."
    conda create -n multiscale2 python=3.9 -y
    echo "请运行: conda activate multiscale2"
    exit 0
fi

echo "[3/6] 安装 PyTorch 生态 (CUDA 11.7)..."
# 注意: dgl 需要从 DGL wheel 索引安装，且使用 dgl 1.1.3 兼容 torch 1.13.x
pip install \
    torch==1.13.1 \
    torchvision==0.14.1 \
    torchaudio==0.13.1 \
    -f https://data.dgl.ai/wheels/cu117/released.html \
    --index-url https://download.pytorch.org/whl/cu117

# dgl 1.1.3 从 DGL wheel 索引安装 (独立安装以避免依赖冲突)
pip install dgl==1.1.3 -f https://data.dgl.ai/wheels/cu117/released.html

pip install \
    e3nn==0.5.1 \
    ml-collections==0.1.1

# se3-transformer: 必须从 huhlim 的 GitHub 仓库安装
pip install git+https://github.com/huhlim/SE3Transformer.git

echo "[4/6] 安装 Conda 包..."
mamba install -c conda-forge -y \
    openmm=8.2.0 \
    numpy=1.24.0 \
    pandas=2.1.1 \
    scipy=1.13.0 \
    matplotlib \
    networkx \
    gromacswrapper \
    parmed 

pip install git+https://github.com/feiglab/mdsim.git


echo "[5/6] 安装 Python 包..."
pip install \
    mdanalysis==2.6.1 \
    biopython==1.81 \
    numba==0.60.0 \
    tqdm \
    pyyaml \
    jinja2 \
    localcider\
    statsmodels

echo "[6/6] 安装 mdtraj..."
pip install mdtraj==1.10.0

# 可选: 从源码安装 mdtraj (huhlim fork)
# pip install git+https://github.com/huhlim/mdtraj.git

echo ""
echo "==========================================="
echo "✅ 安装完成!"
echo "==========================================="
echo ""
echo "验证安装:"
echo "  python -c \"import torch; print(f'PyTorch: {torch.__version__}, CUDA: {torch.version.cuda}')\""
echo "  python -c \"import openmm; print(f'OpenMM: {openmm.__version__}')\""
echo "  python -c \"import MDAnalysis; print(f'MDAnalysis: {MDAnalysis.__version__}')\""
echo ""
echo "如果遇到问题，请运行:"
echo "  python scripts/verify_env.py"

