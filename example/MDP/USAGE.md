# MDP任务使用指南

## 概述

MDP (Membrane-Disrupting Peptides) 任务用于模拟蛋白质与膜的相互作用或在特定区域内具有明确结构的蛋白质。这个示例演示如何使用 `multiscale2` 软件运行CALVADOS粗粒化模拟。

## 文件说明

- `config.yaml` - 主配置文件，包含所有模拟参数
- `TDP43.pdb` - 蛋白质结构文件 (可替换为自己的PDB)
- `domains.yaml` - 蛋白质结构域定义文件，用于施加内部距离限制
- `run_example.py` - 交互式运行脚本
- `USAGE.md` - 本使用说明文档

**注意**: `components.yaml` 文件将由CALVADOS的 `prepare.py` 脚本自动生成，不需要预先提供。

## 快速开始

### 1. 基本运行
```bash
cd example/MDP
conda activate cg2all
export CUDA_VISIBLE_DEVICES=0
python run_example.py
```

## 配置说明

### config.yaml 主要参数

```yaml
# 输入文件配置
input_files:
  structure_pdb: "TDP43.pdb"  # 蛋白质结构文件

# 蛋白质配置
protein:
  name: "TDP43"         # PDB文件名中的蛋白质名称 (不含.pdb)
  nmol: 10              # 蛋白质分子数量

# CALVADOS模拟参数
cg_calvados:
  task_type: "MDP"
  fdomains: "domains.yaml"    # 结构域定义文件
  slab_width: 40              # slab拓扑宽度 (nm)
  steps: 50000000             # 模拟步数
  box: [20.0, 20.0, 80.0]     # 模拟盒子大小 (nm)
```

### 自定义配置

#### 替换蛋白质结构
1.  将你的蛋白质结构保存为PDB格式 (例如, `MyProtein.pdb`) 并放入 `example/MDP` 目录。
2.  创建一个新的 `domains.yaml` 文件，定义蛋白质的结构域。
3.  在 `config.yaml` 中更新以下字段：
    ```yaml
    input_files:
      structure_pdb: "MyProtein.pdb"
    
    protein:
      name: "MyProtein"
      nmol: 20  # 根据需要调整数量
    
    cg_calvados:
      fdomains: "my_domains.yaml"
    ```

## 工作流程

1.  **准备阶段**: 生成器读取 `config.yaml`，找到指定的PDB和domains文件。
2.  **环境设置**: 将PDB, domains和CALVADOS残基定义文件复制到输出目录的 `input` 子目录中。
3.  **脚本生成**: 生成 `prepare.py` 脚本，该脚本引用 `input` 目录中的文件。
4.  **模拟阶段**: `prepare.py` 和 `run.py` 脚本执行模拟。

## 故障排除

### 常见问题

1.  **文件未找到 (FileNotFoundError)**: 确保 `config.yaml` 中指定的 `structure_pdb` 和 `fdomains` 文件名正确无误，并且这些文件与 `config.yaml` 位于同一目录中。
2.  **CALVADOS错误**: 检查 `domains.yaml` 中定义的残基范围是否与PDB文件中的实际残基数量匹配。
