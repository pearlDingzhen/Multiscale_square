# IDR任务使用指南

## 概述

IDR (Intrinsically Disordered Regions) 任务用于模拟内在无序蛋白质区域的行为。这个示例演示如何使用 `multiscale2` 软件运行CALVADOS粗粒化模拟。

## 文件说明

- `config.yaml` - 主配置文件，包含所有模拟参数
- `protein.fasta` - 蛋白质序列文件（可替换为自己的序列）
- `run_example.py` - 交互式运行脚本
- `USAGE.md` - 本使用说明文档

**注意**: `components.yaml` 文件将由CALVADOS的 `prepare.py` 脚本自动生成，不需要预先提供。

## 快速开始

### 1. 基本运行
```bash
cd example/IDR
conda activate cg2all
export CUDA_VISIBLE_DEVICES=0
python run_example.py
```

### 2. 手动运行
```bash
cd example/IDR
python -m multiscale2.calvados_generator config.yaml --output_dir my_idr_simulation --gpu_id 0 --replica 1
```

## 配置说明

### config.yaml 主要参数

```yaml
# 蛋白质配置
protein:
  name: "FUS_LC"        # 从FASTA文件中选择的序列名称
  nmol: 20              # 蛋白质分子数量

# CALVADOS模拟参数
cg_calvados:
  task_type: "IDR"           # 任务类型
  steps: 100000             # 模拟步数
  wfreq: 1000               # 轨迹输出频率
  temp: 310.0               # 温度 (K)
  ionic: 0.15               # 离子浓度 (M)
  box: [25.0, 25.0, 30.0]   # 模拟盒子大小 (nm)
  topol: "slab"             # 拓扑类型
```

### 自定义配置

#### 修改模拟参数
```yaml
cg_calvados:
  steps: 500000      # 增加模拟步数
  wfreq: 5000        # 调整输出频率
  temp: 300.0        # 修改温度
  ionic: 0.1         # 修改离子浓度
  box: [30.0, 30.0, 40.0]  # 增大盒子
```

#### 修改蛋白质配置
```yaml
protein:
  name: "TDP43_FULL"     # 选择不同的蛋白质序列
  nmol: 50               # 增加蛋白质数量
```

#### 替换蛋白质序列
1. 将你的蛋白质序列保存为FASTA格式
2. 替换 `protein.fasta` 文件
3. 在 `config.yaml` 中更新 `protein.name` 为FASTA文件中的序列名称
4. 调整 `protein.nmol` 为所需的蛋白质数量

**示例**: 如果FASTA文件包含多个序列：
```fasta
>FUS_LC
MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS
>TDP43_FULL
MNIDVLEKNVTLSKKDKEDFLIWLLSAGSTFSTSLGSGMMGMLGSSLTDLQSFLIDASLQAETCQLQFTLEKVCQAQGKLLDQCWLQQDQDAKIGIFKGVGACVKIAQDKNTVEKQGQSKMSGFLLFKDPTTMTCGTVEQQNQIRLVEGNLGAYTSSLSSIRASLTFNRGFVAGCNVETIRWFLRGELLDSCVNIFHLQDVFSEHSGQALVLSHGFKITPELVSSLTPSLPTLDPFLTDITPFDGIDPFRLLQDGCWTADLVVHNNLKDLSNVQKPRPFAAGLHEGYCVHSIAQTNLRIPLRRIPFV
```

在 `config.yaml` 中可以选择：
```yaml
protein:
  name: "FUS_LC"        # 或 "TDP43_FULL"
  nmol: 20
```

## 输出文件

运行完成后，会在输出目录中生成：

### 主要输出文件
- `*.dcd` - 轨迹文件 (主要输出，包含所有模拟帧)
- `*.pdb` - 结构文件
- `restart.chk` - 重启检查点文件

### 脚本文件
- `prepare.py` - CALVADOS准备脚本
- `run.py` - CALVADOS运行脚本
- `config.yaml` - CALVADOS配置文件
- `components.yaml` - CALVADOS组件文件 (自动生成)

### 日志文件
- `*.log` - 模拟日志文件

## 工作流程

1. **准备阶段**: 生成器读取配置文件，创建CALVADOS脚本
2. **生成阶段**: `prepare.py` 自动生成 `components.yaml` 和其他必要文件
3. **模拟阶段**: `run.py` 执行CALVADOS模拟
4. **输出阶段**: 生成轨迹文件和分析结果

## 性能优化建议

### 系统大小调整
- **小系统** (< 1000个粒子): 可以使用默认参数
- **中等系统** (1000-5000个粒子): 建议增加模拟步数到500k-1M
- **大系统** (> 5000个粒子): 考虑使用多副本并行

### GPU内存优化
- 根据GPU内存调整蛋白质数量
- 大系统可能需要减少 `wfreq` 来节省内存

### 模拟时间优化
- 使用 `restart.chk` 文件进行续算
- 考虑使用多副本并行运行

## 故障排除

### 常见问题

1. **CUDA错误**: 确保CUDA环境正确设置
2. **内存不足**: 减少蛋白质数量或增加 `wfreq`
3. **轨迹文件为空**: 检查模拟是否正常完成
4. **路径错误**: 确保使用相对路径

### 调试步骤

1. 检查环境变量: `echo $CUDA_VISIBLE_DEVICES`
2. 验证CALVADOS安装: `python -c "import calvados"`
3. 检查GPU可用性: `nvidia-smi`
4. 查看日志文件了解详细错误信息

## 下一步

完成IDR模拟后，可以：

1. **分析轨迹**: 使用CALVADOS分析工具
2. **可视化**: 使用VMD或PyMOL查看结构
3. **继续工作流**: 进行反向映射和全原子模拟
4. **参数优化**: 根据结果调整模拟参数
