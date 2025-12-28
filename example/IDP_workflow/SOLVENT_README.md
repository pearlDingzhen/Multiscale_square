# Solvent Preparation Workflow

这个脚本实现了IDP系统的显式溶剂体系构建和优化工作流，包含三个主要步骤：

## 功能概述

### 三步工作流

1. **真空蛋白质优化** - 在真空中优化蛋白质结构
2. **构建显式溶剂体系** - 添加溶剂分子和离子（可选抗冻水）
3. **显式溶剂盒子优化** - 优化完整的溶剂化体系

### 关键特性

- **自动化的三步流程**：从OpenMM输出到完整的溶剂化体系
- **智能参数设置**：根据config.yaml自动配置所有参数
- **抗冻水支持**：可选择性添加抗冻水分子
- **离子浓度控制**：自动添加指定浓度的离子
- **不同的emtol设置**：
  - 真空优化：steep=500, cg=100
  - 溶剂优化：steep=2000, cg=200

## 使用方法

### 基本用法

```bash
# 在example/IDP_workflow目录中运行
python run_solvent.py [config.yaml] [openmm_output_dir]
```

### 示例

```bash
# 在example/IDP_workflow目录中运行
cd example/IDP_workflow

# 使用默认配置（推荐）
python run_solvent.py

# 或者指定参数
python run_solvent.py config.yaml output_openmm/openmm
```

### 参数说明

- `config.yaml`: 包含所有模拟参数的配置文件
- `openmm_output_dir`: OpenMM优化输出的目录路径（包含conf.gro和PACE.top）

## 配置文件要求

脚本会从config.yaml中读取以下参数：

```yaml
protein:
  name: "FUS_LC"        # 蛋白质名称
  nmol: 20              # 分子数量

gromacs_explicit_solvent_build_and_equilibration:
  num_cpus: 24                    # CPU数量
  ions_concentration: 0.15        # 离子浓度 (M)
  temperature: 300.0              # 温度 (K)
  pressure: 1.0                   # 压力 (bar)
  Use_anti_freeze_water: True     # 是否使用抗冻水
  Anti_free_water_ratio: 0.1      # 抗冻水比例
```

## 输出文件结构

```
output_solvent/
├── vacuum_optimization/          # 真空优化
│   ├── conf.gro                 # 输入结构
│   ├── PACE.top                 # 拓扑文件
│   ├── em_steep.mdp             # Steep descent参数
│   ├── em_cg.mdp                # Conjugate gradient参数
│   ├── em_steep.gro             # Steep descent输出
│   └── em_cg.gro                # 最终真空优化结构
└── explicit_solvent/             # 显式溶剂
    ├── system.gro               # 溶剂化结构
    ├── topol.top                # 溶剂化拓扑
    ├── em_steep.mdp             # 溶剂steep descent参数
    ├── em_cg.mdp                # 溶剂conjugate gradient参数
    ├── em_steep.gro             # 溶剂steep descent输出
    └── em_cg.gro                # 最终优化结构
```

## 工作流步骤详解

### Step 1: 真空蛋白质优化

- 从OpenMM输出获取`conf.gro`和`PACE.top`
- 在真空中进行能量极小化
- 使用emtol=500 (steep) 和 emtol=100 (cg)
- 输出：`vacuum_optimization/em_cg.gro`

### Step 2: 构建显式溶剂体系

- 使用GROMACS的`solvate`命令添加溶剂
- 可选添加抗冻水（ASOL分子）
- 添加指定浓度的离子
- 输出：`explicit_solvent/system.gro`

### Step 3: 显式溶剂盒子优化

- 对完整的溶剂化体系进行能量极小化
- 使用emtol=2000 (steep) 和 emtol=200 (cg)
- 输出：`explicit_solvent/em_cg.gro`

## 依赖要求

- GROMACS (gmx命令)
- Python 3.7+
- 依赖包：yaml, pathlib, shutil, subprocess
- multiscale2包中的utils模块

## 注意事项

1. **文件路径**：确保OpenMM输出目录包含`conf.gro`和`PACE.top`文件
2. **GROMACS环境**：确保GROMACS已正确安装并配置环境变量
3. **计算资源**：根据系统大小调整CPU数量
4. **存储空间**：溶剂化体系会产生大量文件，确保有足够存储空间

## 故障排除

### 常见问题

1. **文件未找到**：检查OpenMM输出目录路径
2. **GROMACS命令失败**：检查GROMACS安装和环境变量
3. **内存不足**：减少CPU数量或增加系统内存
4. **拓扑文件错误**：检查PACE.top文件格式

### 调试模式

脚本会输出详细的执行日志，包括：
- 文件复制操作
- GROMACS命令执行
- 错误信息和堆栈跟踪

## 后续步骤

完成溶剂准备后，系统可以用于：
1. 平衡模拟（NVT/NPT）
2. 生产模拟
3. 分析计算

使用`multiscale2.gromacs_equilibrate.IDPGromacsEquilibrator`进行后续的平衡和生产模拟。
