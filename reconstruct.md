一、 Multiscale² 重构蓝图 (Refactoring Blueprint)
重构将分为四个阶段，旨在实现从“生成脚本”到“执行逻辑”的范式转变。

1. 统一命令中心 (Unified CLI Engine)
目标：模仿 GROMACS，提供 ms2 <subcommand> 接口。

逻辑：放弃 WorkflowGenerator 生成 Python 脚本的做法。改为由 CLI 直接调用内部封装好的 Manager 类。

价值：用户不再需要面对一堆生成的 .py 文件，只需维护一个 config.yaml。

2. 状态化项目管理 (Stateful Project Context)
目标：引入 ms2_state.json 记录项目进度。

逻辑：创建一个 ProjectContext 类，它记录了当前工作的目录、已完成的阶段、每一步生成的关键文件路径（如 checkpoint.pdb 或 PACE.top）。

价值：实现断点续传。如果全原子平衡崩溃，用户修复配置后只需重新执行 ms2 aa，系统会自动找到 Stage 3 生成的 PACE 结构。

3. 静态依赖集成 (Static Vendoring & Isolation)
目标：彻底解决 cg2all 与 mdtraj 的依赖冲突。

逻辑：将 cg2all 和 CALVADOS 的核心算法代码放入 extern/ 目录。

深度优化：重写 cg2all 中依赖魔改 mdtraj 的 IO 部分，统一改用 MDAnalysis 进行原子坐标读写。

4. 物理缓冲层自动化 (Automated PACE "Soft Landing")
目标：强化 PACE 在全原子转换中的缓冲作用。

逻辑：在 PaceManager 中实现自动化的“渐进式约束”：

Position Restraints：初次 EM 时强约束所有重原子。

Sidechain Release：释放侧链进行 PACE 优化。

Backbone Release：最后释放骨架，确保结构在进入全原子硬势能面前已经物理合理。

二、 软件核心目录结构

```
Multiscale_square/              # 项目根目录
├── multiscale2/                # 主程序包
│   ├── __init__.py
│   ├── cli/                    # 命令行接口
│   │   ├── main.py             # 入口
│   │   └── commands/           # 子命令
│   ├── core/                   # 核心引擎
│   │   ├── context.py          # 项目状态管理
│   │   ├── config.py           # 配置校验
│   │   └── executor.py         # 子进程封装
│   ├── stages/                 # 模拟阶段
│   │   ├── base.py
│   │   ├── cg_stage.py         # CALVADOS
│   │   ├── backmap_stage.py    # cg2all
│   │   ├── pace_stage.py       # PACE
│   │   └── aa_stage.py         # 全原子平衡
│   ├── extern/                 # 静态集成第三方包 ⬅️ 重点
│   │   ├── __init__.py
│   │   ├── ms2_cg2all/         # cg2all 全原子重建
│   │   │   ├── __init__.py     # 暴露 convert_cg2all
│   │   │   ├── lib/            # 核心算法
│   │   │   │   ├── libpdb.py   # PDB 处理（已修改：移除 bfactors）
│   │   │   │   ├── libcg.py    # CG 模型（已修改：移除 bfactors）
│   │   │   │   ├── libdata.py
│   │   │   │   ├── libmodel.py
│   │   │   │   ├── libloss.py
│   │   │   │   └── snippets.py
│   │   │   ├── model/          # 模型权重
│   │   │   │   ├── CalphaBasedModel.ckpt
│   │   │   │   └── ...
│   │   │   └── extern/         # cg2all 的依赖
│   │   │       └── ms2_se3transformer/
│   │   ├── ms2_calvados/       # CALVADOS 粗粒化 ⬅️ 待集成
│   │   │   ├── __init__.py     # 暴露子模块
│   │   │   ├── calvados/       # ⬅️ 原始 CALVADOS 源码（复制自 source_code_of_extern/CALVADOS-main/）
│   │   │   │   ├── __init__.py # 导入 cfg, sim, build, interactions 等
│   │   │   │   ├── cfg.py      # Config, Job, Components 类
│   │   │   │   ├── sim.py      # 模拟引擎（OpenMM）
│   │   │   │   ├── build.py    # 系统构建
│   │   │   │   ├── interactions.py  # 相互作用势
│   │   │   │   ├── analysis.py # 分析工具
│   │   │   │   ├── sequence.py # 序列处理
│   │   │   │   ├── utilities.py
│   │   │   │   ├── components.py
│   │   │   │   └── data/       # 残基参数
│   │   │   │       ├── residues.csv
│   │   │   │       └── templates/
│   │   │   └── wrapper.py      # 包装器（保持 CalvadosGenerator 接口）
│   │   ├── ms2_cocomo/         # 可选：cocomo
│   │   ├── ms2_Mpipi-re/       # 可选：mpipi-recharged
│   │   └── ms2_openabc/        # 可选：openabc (HPS/MOFF)
│   ├── ff/                     # 力场参数
│   │   └── pace/               # PACE 参数文件
│   ├── assets/                 # 静态资源
│   │   ├── templates/          # GROMACS mdp 模板
│   │   └── residues/           # 残基定义
│   └── utils/                  # 通用工具
├── tests/                      # 测试
├── examples/                   # 示例
├── docs/                       # 文档
├── requirements.txt
└── setup.py
```

### extern/ 目录说明

这是我们「静态集成」策略的核心，所有第三方算法源码直接嵌入，不修改上游代码，仅做兼容性适配。

| 包名 | 功能 | 修改状态 |
|------|------|----------|
| ms2_cg2all | 全原子重建 | ✅ 已移除 bfactors 依赖 |
| ms2_calvados | 粗粒化模拟 | 🚧 待集成（下一步） |
| ms2_cocomo | 可选：粗粒化 | ⏳ 待定 |
| ms2_Mpipi-re | 可选：π-π 相互作用 | ⏳ 待定 |
| ms2_openabc | 可选：HPS/MOFF | ⏳ 待定 |

---

### ms2_calvados 集成计划

**来源**：`source_code_of_extern/CALVADOS-main/`

**核心接口**（需要保持兼容）：
```python
from calvados.cfg import Config, Job, Components
```

**集成步骤**：

1. **复制源码**
   ```bash
   cp -r source_code_of_extern/CALVADOS-main/calvados \
         multiscale2/extern/ms2_calvados/
   ```

2. **创建包装器** (`multiscale2/extern/ms2_calvados/wrapper.py`)
   ```python
   # 保持 CalvadosGenerator 接口
   from .calvados.cfg import Config, Components
   
   class CalvadosGenerator:
       def __init__(self, config_path): ...
       def generate_and_run(self, output_dir, protein_name, gpu_id, replica): ...
   ```

3. **更新导入**
   - 移除 `from calvados.cfg` → 改为 `from ms2_calvados.calvados.cfg`
   - 保持 `residues.csv` 在 `data/` 目录

**需要修改的部分**：
- setup.py 中的依赖声明（改为静态集成后不需要）
- 可能需要调整 data/ 目录的路径引用

### 迁移计划

**已完成**：
- multiscale2/extern/ms2_cg2all/lib/libpdb.py - 移除 bfactors
- multiscale2/extern/ms2_cg2all/lib/libcg.py - 移除 bfactors

**下一步**：
1. 迁移 multiscale2/*.py → multiscale2/core/
2. 迁移 workflow_stages/*.py → multiscale2/stages/
3. 集成 ms2_calvados