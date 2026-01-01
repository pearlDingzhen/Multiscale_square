# CG Simulation 模块代码审查报告

**文件**: `multiscale2/src/cg.py`  
**审查日期**: 2024年12月30日  
**代码行数**: 923行

---

## 一、大体分析

### 1.1 文件定位与核心功能

`cg.py` 是 multiscale2 项目的核心模块，提供统一的粗粒化（Coarse-Grained, CG）分子动力学模拟接口。该模块设计为支持多种力场，包括 CALVADOS、HPS、MOFF、COCOMO 和 OpenMpipi。

该文件的核心职责是：
- **配置管理**: 定义模拟配置的数据类（CGSimulationConfig, CGComponent, SimulationParams）
- **模拟器**: 实现统一的模拟器类 CGSimulator
- **多力场支持**: 为不同力场提供独立的 runner 方法
- **结果管理**: 定义模拟结果的数据结构

### 1.2 架构设计

```
cg.py
├── Enums (枚举类)
│   ├── ComponentType: 组件类型（IDP/MDP）
│   ├── TopologyType: 拓扑类型（CUBIC/SLAB）
│   └── Platform: 计算平台（CPU/CUDA）
│
├── Configuration Classes (配置类)
│   ├── SimulationParams: 核心模拟参数
│   ├── CGComponent: 单个组件规格
│   └── CGSimulationConfig: 完整模拟配置
│
├── Simulation Result (结果类)
│   └── SimulationResult: 模拟结果数据类
│
└── CG Simulator (模拟器类)
    └── CGSimulator: 粗粒化模拟器（含多个 runner 方法）
        ├── Setup Methods (设置方法)
        │   ├── __init__()
        │   ├── setup()
        │   ├── prepare_calvados_output()
        │   ├── _copy_input_files()
        │   ├── _ensure_setup()
        │   └── _ensure_not_running()
        │
        ├── Runner Methods (运行方法)
        │   ├── run_calvados()
        │   ├── run_hps()
        │   ├── run_moff()
        │   ├── run_cocomo()
        │   └── run_openmpipi()
        │
        └── Utility Methods (工具方法)
            ├── get_result()
            └── cleanup()
```

### 1.3 代码风格评价

**优点**:
- 使用 dataclass 定义配置类，代码简洁
- 完整的类型注解
- 详细的文档字符串
- 清晰的模块划分
- 良好的枚举设计

**需改进**:
- 部分 runner 方法未完全实现（占位符实现）
- 错误处理可以更细化
- 缺少单元测试覆盖

---

## 二、枚举类详细解读

### 2.1 ComponentType 枚举

**功能**: 定义组件类型枚举。

**枚举值**:
- `IDP = "idp"`: 无序蛋白（Intrinsically Disordered Protein），基于序列模拟
- `MDP = "mdp"`: 折叠蛋白（Molecular Dynamics Protein），基于结构模拟

**代码质量**: ⭐⭐⭐⭐⭐
- 简洁明了
- 描述清晰

**代码示例**:
```python
class ComponentType(Enum):
    """组件类型"""
    IDP = "idp"   # 无序蛋白 - 基于序列
    MDP = "mdp"   # 折叠蛋白 - 基于结构
```

### 2.2 TopologyType 枚举

**功能**: 定义系统拓扑类型枚举。

**枚举值**:
- `CUBIC = "cubic"`: 立方体盒子，使用网格放置
- `SLAB = "slab"`: 平面限制，用于相分离模拟

**代码质量**: ⭐⭐⭐⭐⭐
- 简洁明了
- 注释清晰

**代码示例**:
```python
class TopologyType(Enum):
    """拓扑类型"""
    CUBIC = "cubic"   # 立方体盒子（网格放置）
    SLAB = "slab"     # 平面限制（相分离）
```

### 2.3 Platform 枚举

**功能**: 定义计算平台枚举。

**枚举值**:
- `CPU = "CPU"`: CPU 计算
- `CUDA = "CUDA"`: GPU 计算

**代码质量**: ⭐⭐⭐⭐
- 简洁明了

**代码示例**:
```python
class Platform(Enum):
    """计算平台"""
    CPU = "CPU"
    CUDA = "CUDA"
```

---

## 三、配置类详细解读

### 3.1 SimulationParams 类

**功能**: 核心模拟参数数据类。

**属性**:
- `_DT`: 时间步长，默认 0.01 ps（10 fs）
- `_FRICTION`: 摩擦系数，默认 0.01
- `steps`: 总积分步数，默认 100000
- `wfreq`: 写入频率，默认 1000
- `platform`: 计算平台，默认 CUDA
- `verbose`: 详细输出，默认 True

**设计特点**:
- dt 和 friction 使用固定默认值，不让用户输入
- 使用 dataclass 的 field 机制定义默认值

**代码质量**: ⭐⭐⭐⭐
- 清晰的默认值设计
- 有 to_dict() 和 from_dict() 方法

**代码示例**:
```python
@dataclass
class SimulationParams:
    """核心模拟参数
    
    Notes:
        dt 和 friction 使用固定默认值，不让用户输入：
        - dt = 0.01 ps (10 fs) - 所有力场的通用默认值
        - friction = 0.01 - OpenMM LangevinMiddleIntegrator 默认值
    """
    # 固定默认值（不让用户输入）
    _DT: float = 0.01       # 时间步长（ps）- 所有力场通用 10fs
    _FRICTION: float = 0.01 # 摩擦系数 - OpenMM 默认值
    
    steps: int = 100000          # 总积分步数
    wfreq: int = 1000            # 写入频率（每N步保存一次）
    platform: Platform = Platform.CUDA
    verbose: bool = True
    
    def to_dict(self) -> Dict:
        d = {
            'steps': self.steps,
            'wfreq': self.wfreq,
            'platform': self.platform.value,
            'verbose': self.verbose,
        }
        return {k: v for k, v in d.items() if v is not None}
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'SimulationParams':
        if 'platform' in d and isinstance(d['platform'], str):
            d['platform'] = Platform(d['platform'])
        d = {k: v for k, v in d.items() if v is not None}
        return cls(**d)
```

**潜在问题**:
1. `_DT` 和 `_FRICTION` 使用下划线前缀，但实际是公共属性
2. to_dict() 中移除了 None 值，但缺少 `_DT` 和 `_FRICTION`

**改进建议**:
```python
@dataclass
class SimulationParams:
    """核心模拟参数"""
    
    # 固定默认值
    DT: float = 0.01       # 时间步长（ps）
    FRICTION: float = 0.01 # 摩擦系数
    
    steps: int = 100000
    wfreq: int = 1000
    platform: Platform = Platform.CUDA
    verbose: bool = True
    
    def to_dict(self) -> Dict:
        """转换为字典"""
        return {
            'dt': self.DT,
            'friction': self.FRICTION,
            'steps': self.steps,
            'wfreq': self.wfreq,
            'platform': self.platform.value if isinstance(self.platform, Platform) else self.platform,
            'verbose': self.verbose,
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'SimulationParams':
        """从字典创建"""
        platform = d.get('platform', 'CUDA')
        if isinstance(platform, str):
            platform = Platform(platform)
        
        return cls(
            DT=d.get('dt', 0.01),
            FRICTION=d.get('friction', 0.01),
            steps=d.get('steps', 100000),
            wfreq=d.get('wfreq', 1000),
            platform=platform,
            verbose=d.get('verbose', True),
        )
```

### 3.2 CGComponent 类

**功能**: 单个组件规格数据类。

**属性**:
- **基础属性**: name, type, nmol
- **输入文件**: ffasta, fpdb, fdomains, fpae
- **约束设置**: restraint, restraint_type, use_com, k_harmonic, colabfold
- **电荷设置**: charge_termini
- **派生属性**: seq, nres

**代码质量**: ⭐⭐⭐⭐
- 属性丰富
- 有完整的序列化和验证方法

**关键方法**:

#### 3.2.1 `to_dict(self) -> Dict`

**功能**: 将组件转换为字典。

**实现逻辑**:
1. 构建基础字典
2. 添加可选文件路径
3. 移除 None 值

**代码示例**:
```python
def to_dict(self) -> Dict:
    d = {
        'name': self.name,
        'type': self.type.value,
        'nmol': self.nmol,
        'restraint': self.restraint,
        'restraint_type': self.restraint_type,
        'use_com': self.use_com,
        'k_harmonic': self.k_harmonic,
        'colabfold': self.colabfold,
        'charge_termini': self.charge_termini,
    }
    if self.ffasta:
        d['ffasta'] = self.ffasta
    if self.fpdb:
        d['fpdb'] = self.fpdb
    if self.fdomains:
        d['fdomains'] = self.fdomains
    if self.fpae:
        d['fpae'] = self.fpae
    return d
```

#### 3.2.2 `from_dict(cls, d: Dict) -> 'CGComponent'`

**功能**: 从字典创建组件。

**实现逻辑**:
1. 解析类型（支持字符串转枚举）
2. 构建组件实例

**代码示例**:
```python
@classmethod
def from_dict(cls, d: Dict) -> 'CGComponent':
    comp_type = d.get('type', 'idp')
    if isinstance(comp_type, str):
        comp_type = ComponentType(comp_type.lower())
    
    return cls(
        name=d['name'],
        type=comp_type,
        nmol=d.get('nmol', 1),
        ffasta=d.get('ffasta'),
        fpdb=d.get('fpdb'),
        fdomains=d.get('fdomains'),
        fpae=d.get('fpae'),
        restraint=d.get('restraint', False),
        restraint_type=d.get('restraint_type', 'harmonic'),
        use_com=d.get('use_com', False),
        k_harmonic=d.get('k_harmonic', 700.0),
        colabfold=d.get('colabfold', 1),
        charge_termini=d.get('charge_termini', 'both'),
    )
```

#### 3.2.3 `validate(self) -> List[str]`

**功能**: 验证组件配置。

**验证逻辑**:
- IDP: 检查 ffasta 文件存在性
- MDP: 检查 fpdb 文件存在性
- 约束: 检查 fdomains 文件存在性（仅文件路径，不检查内联 YAML）

**代码质量**: ⭐⭐⭐⭐⭐
- 健壮的验证逻辑
- 支持内联 YAML 格式

**代码示例**:
```python
def validate(self) -> List[str]:
    """验证配置
    
    fdomains 支持两种格式：
    1. 文件路径：'domains.yaml' - 检查文件是否存在
    2. 内联 YAML：'TDP43:\n  - [3, 76]\n...' - 不检查
    """
    errors = []
    
    def _is_inline_yaml(text: str) -> bool:
        """检查是否是内联 YAML"""
        if not text:
            return False
        stripped = text.strip()
        if stripped.startswith('{') or stripped.startswith('['):
            return True
        if '\n' in stripped and (':' in stripped or stripped.startswith('-')):
            return True
        return False
    
    if self.type == ComponentType.IDP:
        if not self.ffasta:
            errors.append(f"Component '{self.name}': IDP requires ffasta file")
        elif not os.path.exists(self.ffasta):
            errors.append(f"Component '{self.name}': FASTA file not found: {self.ffasta}")
    elif self.type == ComponentType.MDP:
        if not self.fpdb:
            errors.append(f"Component '{self.name}': MDP requires fpdb file")
        elif not os.path.exists(self.fpdb):
            errors.append(f"Component '{self.name}': PDB file not found: {self.fpdb}")
        if self.restraint and self.fdomains:
            if not _is_inline_yaml(self.fdomains) and not os.path.exists(self.fdomains):
                errors.append(f"Component '{self.name}': Domains file not found: {self.fdomains}")
    return errors
```

### 3.3 CGSimulationConfig 类

**功能**: 完整模拟配置数据类。

**属性**:
- **系统信息**: system_name
- **环境**: box, temperature, ionic
- **拓扑**: topol
- **模拟参数**: simulation
- **组件列表**: components
- **输出**: output_dir
- **元数据**: config_path, created_at

**代码质量**: ⭐⭐⭐⭐
- 完整的配置结构
- 丰富的辅助方法

**关键方法**:

#### 3.3.1 `__init__` 和属性设置

**功能**: 初始化配置。

**默认值**:
- system_name: "cg_simulation"
- box: [25.0, 25.0, 30.0]
- temperature: 310.0 K
- ionic: 0.15 M
- topol: CUBIC
- simulation: SimulationParams()
- output_dir: "output_cg"

**代码示例**:
```python
@dataclass
class CGSimulationConfig:
    """完整模拟配置
    
    示例 YAML 结构：
        system_name: my_simulation
        box: [25.0, 25.0, 30.0]
        temperature: 310.0
        ionic: 0.15
        topol: cubic  # 或 slab
        
        simulation:
          steps: 100000
          wfreq: 1000
          platform: CUDA
        
        components:
          - name: protein_A
            type: IDP
            nmol: 20
            ffasta: input/protein_A.fasta
    """
    # 系统信息
    system_name: str = "cg_simulation"
    
    # 环境
    box: List[float] = field(default_factory=lambda: [25.0, 25.0, 30.0])
    temperature: float = 310.0       # Kelvin
    ionic: float = 0.15              # Molar
    
    # 拓扑
    topol: TopologyType = TopologyType.CUBIC
    
    # 模拟参数
    simulation: SimulationParams = field(default_factory=SimulationParams)
    
    # 组件列表
    components: List[CGComponent] = field(default_factory=list)
    
    # 输出
    output_dir: str = "output_cg"
    
    # 元数据
    config_path: Optional[str] = None
    created_at: str = field(default_factory=lambda: str(__import__('datetime').datetime.now()))
```

#### 3.3.2 `add_component(self, component: CGComponent)`

**功能**: 添加组件。

**代码示例**:
```python
def add_component(self, component: CGComponent):
    """添加组件"""
    self.components.append(component)
```

#### 3.3.3 `get_component(self, name: str) -> Optional[CGComponent]`

**功能**: 根据名称获取组件。

**代码示例**:
```python
def get_component(self, name: str) -> Optional[CGComponent]:
    """根据名称获取组件"""
    for comp in self.components:
        if comp.name == name:
            return comp
    return None
```

#### 3.3.4 `total_molecules(self) -> int`

**功能**: 计算总分子数。

**代码示例**:
```python
def total_molecules(self) -> int:
    """计算总分子数"""
    return sum(comp.nmol for comp in self.components)
```

#### 3.3.5 `validate(self) -> List[str]`

**功能**: 验证完整配置。

**验证内容**:
- system_name 非空
- box 为 3 元素列表
- 至少一个组件
- 所有组件验证通过

**代码示例**:
```python
def validate(self) -> List[str]:
    """验证配置"""
    errors = []
    if not self.system_name:
        errors.append("system_name is required")
    if len(self.box) != 3:
        errors.append("box must be a list of 3 values [x, y, z]")
    if not self.components:
        errors.append("At least one component is required")
    for comp in self.components:
        errors.extend(comp.validate())
    return errors
```

#### 3.3.6 `to_dict(self) -> Dict`

**功能**: 转换为字典。

**代码示例**:
```python
def to_dict(self) -> Dict:
    """转换为字典"""
    return {
        'system_name': self.system_name,
        'box': self.box,
        'temperature': self.temperature,
        'ionic': self.ionic,
        'topol': self.topol.value if isinstance(self.topol, TopologyType) else self.topol,
        'simulation': self.simulation.to_dict(),
        'components': [c.to_dict() for c in self.components],
        'output_dir': self.output_dir,
    }
```

#### 3.3.7 `to_yaml(self, path: str = None)`

**功能**: 保存到 YAML 文件。

**代码示例**:
```python
def to_yaml(self, path: str = None):
    """保存到 YAML 文件"""
    d = self.to_dict()
    if path:
        with open(path, 'w') as f:
            yaml.dump(d, f, default_flow_style=False, sort_keys=False)
    else:
        return yaml.dump(d, default_flow_style=False, sort_keys=False)
```

#### 3.3.8 `from_dict(cls, d: Dict) -> 'CGSimulationConfig'`

**功能**: 从字典创建配置。

**代码示例**:
```python
@classmethod
def from_dict(cls, d: Dict) -> 'CGSimulationConfig':
    """从字典创建"""
    topol = d.get('topol', 'cubic')
    if isinstance(topol, str):
        topol = TopologyType(topol.lower())
    
    sim_dict = d.get('simulation', {})
    if isinstance(sim_dict, dict):
        simulation = SimulationParams.from_dict(sim_dict)
    else:
        simulation = SimulationParams()
    
    components = []
    for comp_dict in d.get('components', []):
        components.append(CGComponent.from_dict(comp_dict))
    
    return cls(
        system_name=d.get('system_name', 'cg_simulation'),
        box=d.get('box', [25.0, 25.0, 30.0]),
        temperature=d.get('temperature', 310.0),
        ionic=d.get('ionic', 0.15),
        topol=topol,
        simulation=simulation,
        components=components,
        output_dir=d.get('output_dir', 'output_cg'),
        config_path=d.get('config_path'),
    )
```

#### 3.3.9 `from_yaml(cls, path: str) -> 'CGSimulationConfig'`

**功能**: 从 YAML 文件加载配置。

**代码示例**:
```python
@classmethod
def from_yaml(cls, path: str) -> 'CGSimulationConfig':
    """从 YAML 文件加载"""
    with open(path, 'r') as f:
        d = yaml.safe_load(f)
    d['config_path'] = path
    return cls.from_dict(d)
```

---

## 四、结果类详细解读

### 4.1 SimulationResult 类

**功能**: 模拟结果数据类。

**属性**:
- success: 是否成功
- output_dir: 输出目录
- trajectory: 轨迹文件路径
- structure: 结构文件路径
- checkpoint: 检查点文件路径
- log: 日志文件路径
- metrics: 指标字典
- errors: 错误列表

**代码质量**: ⭐⭐⭐⭐
- 完整的字段定义
- 使用 field 定义默认值

**代码示例**:
```python
@dataclass
class SimulationResult:
    """模拟结果"""
    success: bool = False
    output_dir: str = ""
    trajectory: Optional[str] = None
    structure: Optional[str] = None
    checkpoint: Optional[str] = None
    log: Optional[str] = None
    metrics: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
```

---

## 五、模拟器类详细解读

### 5.1 CGSimulator 类

**功能**: 粗粒化模拟器，提供统一的模拟接口。

**属性**:
- config: 模拟配置
- output_dir: 输出目录
- is_setup: 是否已完成设置
- is_running: 是否正在运行
- _result: 最近的结果

**设计特点**:
- 支持多种力场（CALVADOS, HPS, MOFF, COCOMO, OpenMpipi）
- 统一的设置和运行流程
- 完善的状态管理

#### 5.1.1 `__init__(self, config: CGSimulationConfig)`

**功能**: 初始化模拟器。

**实现逻辑**:
1. 保存配置
2. 初始化状态
3. 验证配置
4. 打印初始化信息

**代码质量**: ⭐⭐⭐⭐
- 清晰的状态初始化
- 验证配置有效性

**代码示例**:
```python
def __init__(self, config: CGSimulationConfig):
    """初始化模拟器
    
    Args:
        config: CGSimulationConfig 实例
    """
    self.config = config
    self.output_dir: Optional[str] = None
    self.is_setup: bool = False
    self.is_running: bool = False
    self._result: Optional[SimulationResult] = None
    
    # 验证配置
    errors = self.config.validate()
    if errors:
        raise ValueError(
            f"Configuration validation failed:\n" +
            "\n".join(f"  - {e}" for e in errors)
        )
    
    print(f"[CGSimulator] Initialized")
    print(f"  System: {config.system_name}")
    print(f"  Components: {len(config.components)}")
    print(f"  Total molecules: {config.total_molecules()}")
```

#### 5.1.2 `setup(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]`

**功能**: 设置模拟环境。

**实现逻辑**:
1. 创建输出目录
2. 复制输入文件
3. 设置 is_setup 标志

**代码质量**: ⭐⭐⭐⭐
- 良好的错误处理
- 支持覆盖模式

**代码示例**:
```python
def setup(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
    """设置模拟环境（通用准备）
    
    创建输出目录并复制输入文件。
    
    Args:
        output_dir: 输出目录
        overwrite: 是否覆盖已存在的目录
    
    Returns:
        生成的文件路径字典
    """
    self._ensure_not_running()
    
    output_dir = os.path.abspath(output_dir)
    
    if os.path.exists(output_dir):
        if not overwrite:
            raise FileExistsError(
                f"Output directory exists: {output_dir}\n"
                f"Use overwrite=True to replace."
            )
    else:
        os.makedirs(output_dir, exist_ok=True)
    
    self.output_dir = output_dir
    
    print(f"\n[CGSimulator] Setting up...")
    print(f"  Output directory: {output_dir}")
    print(f"  System: {self.config.system_name}")
    
    # 复制输入文件
    self._copy_input_files(output_dir)
    
    self.is_setup = True
    print(f"  ✓ Setup complete")
    
    return {
        'output_dir': output_dir,
        'config': os.path.join(output_dir, 'config.yaml'),
    }
```

#### 5.1.3 `prepare_calvados_output(self) -> Dict[str, str]`

**功能**: 准备 CALVADOS 输出的目录结构。

**输出结构**:
```
{output_dir}/
├── {system_name}_CG/
│   ├── raw/                  # 原生输出
│   ├── trajectory.dcd        # 整理后的轨迹
│   ├── final.pdb             # 整理后的最终结构
│   └── simulation.log        # 高层级日志
```

**代码质量**: ⭐⭐⭐⭐
- 清晰的目录结构设计
- 自动备份旧结果

**代码示例**:
```python
def prepare_calvados_output(self) -> Dict[str, str]:
    """准备 CALVADOS 输出的目录结构
    
    统一输出结构：
    {output_dir}/
    ├── {system_name}_CG/
    │   ├── raw/                  # 原生输出
    │   ├── trajectory.dcd        # 整理后的轨迹
    │   ├── final.pdb             # 整理后的最终结构
    │   └── simulation.log        # 高层级日志
    
    Returns:
        包含输出路径的字典
    """
    self._ensure_setup()
    self._ensure_not_running()
    
    expected_suffix = f"{self.config.system_name}_CG"
    if self.output_dir.endswith(expected_suffix):
        output_dir = self.output_dir
        task_name = expected_suffix
    else:
        task_name = expected_suffix
        output_dir = os.path.join(self.output_dir, task_name)
    
    raw_dir = os.path.join(output_dir, 'raw')
    
    # 如果目录已存在，备份后重建
    import shutil
    from datetime import datetime
    
    if os.path.exists(output_dir):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_dir = f"{output_dir}_backup_{timestamp}"
        shutil.move(output_dir, backup_dir)
        print(f"  📁 备份旧结果到: {backup_dir}")
    
    os.makedirs(raw_dir, exist_ok=True)
    
    return {
        'output_dir': output_dir,
        'raw_dir': raw_dir,
        'task_name': task_name,
    }
```

#### 5.1.4 `_copy_input_files(self, output_dir: str)`

**功能**: 复制输入文件到输出目录。

**实现逻辑**:
1. 创建 input 子目录
2. 复制 IDP 的 FASTA 文件
3. 复制 MDP 的 PDB、domains、PAE 文件

**代码质量**: ⭐⭐⭐⭐
- 清晰的输入文件管理
- 使用 shutil.copy2 保持文件元数据

**代码示例**:
```python
def _copy_input_files(self, output_dir: str):
    """复制输入文件到输出目录"""
    input_dir = os.path.join(output_dir, 'input')
    os.makedirs(input_dir, exist_ok=True)
    
    for comp in self.config.components:
        if comp.type == ComponentType.IDP and comp.ffasta:
            if os.path.exists(comp.ffasta):
                shutil.copy2(comp.ffasta, os.path.join(input_dir, os.path.basename(comp.ffasta)))
        
        elif comp.type == ComponentType.MDP:
            if comp.fpdb and os.path.exists(comp.fpdb):
                shutil.copy2(comp.fpdb, os.path.join(input_dir, os.path.basename(comp.fpdb)))
            if comp.fdomains and os.path.exists(comp.fdomains):
                shutil.copy2(comp.fdomains, os.path.join(input_dir, os.path.basename(comp.fdomains)))
            if comp.fpae and os.path.exists(comp.fpae):
                shutil.copy2(comp.fpae, os.path.join(input_dir, os.path.basename(comp.fpae)))
```

#### 5.1.5 `_ensure_setup(self)` 和 `_ensure_not_running(self)`

**功能**: 状态检查方法。

**代码质量**: ⭐⭐⭐⭐
- 清晰的错误消息

**代码示例**:
```python
def _ensure_setup(self):
    """确保已完成设置"""
    if not self.is_setup:
        raise RuntimeError("Simulation not set up. Call setup() first.")

def _ensure_not_running(self):
    """确保未在运行"""
    if self.is_running:
        raise RuntimeError("Simulation is already running")
```

#### 5.1.6 `run_calvados(self, gpu_id: int = 0, **kwargs) -> SimulationResult`

**功能**: 运行 CALVADOS 模拟。

**实现逻辑**:
1. 准备输出目录
2. 调用 CalvadosWrapper 写入配置
3. 运行模拟
4. 整理输出
5. 写入日志

**代码质量**: ⭐⭐⭐⭐
- 清晰的运行流程
- 错误处理完善

**问题**:
1. `elapsed` 变量定义为 0，未实际计时
2. 调用 `wrapper._write_to_dir(raw_dir)` 使用了私有方法

**代码示例**:
```python
def run_calvados(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
    """运行 CALVADOS 模拟
    
    直接委托给 CalvadosWrapper 进行 CALVADOS 模拟。
    
    Args:
        gpu_id: GPU 设备 ID
        **kwargs: 额外参数
    
    Returns:
        SimulationResult
    """
    from .calvados_wrapper import CalvadosWrapper
    import shutil
    import time
    from datetime import datetime
    
    self._ensure_setup()
    self._ensure_not_running()
    
    # 准备输出目录
    dirs = self.prepare_calvados_output()
    output_dir = dirs['output_dir']
    raw_dir = dirs['raw_dir']
    task_name = dirs['task_name']
    
    self.is_running = True
    result = SimulationResult()
    result.output_dir = output_dir
    
    try:
        print(f"\n[CALVADOS] Running simulation via CGSimulator...")
        print(f"  GPU ID: {gpu_id}")
        print(f"  Task: {task_name}")
        print(f"  Raw output: {raw_dir}")
        print(f"  Topology: {self.config.topol.value}")
        
        # 设置 GPU
        os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
        
        # 调用 CalvadosWrapper
        wrapper = CalvadosWrapper(self.config)
        wrapper._write_to_dir(raw_dir)
        
        # 运行模拟
        from multiscale2.extern.ms2_calvados.calvados import sim as calvados_sim
        calvados_sim.run(
            path=raw_dir,
            fconfig='config.yaml',
            fcomponents='components.yaml'
        )
        
        # 组织输出
        self._organize_calvados_output(raw_dir, output_dir, task_name)
        
        # 写入日志
        self._write_simulation_log(output_dir, task_name, 0, True)
        
        result.success = True
        print(f"  ✓ CALVADOS simulation completed")
    
    except Exception as e:
        result.success = False
        result.errors.append(str(e))
        print(f"  ✗ CALVADOS simulation failed: {e}")
    
    finally:
        self.is_running = False
    
    # 设置结果文件路径
    result.trajectory = os.path.join(output_dir, 'trajectory.dcd')
    result.structure = os.path.join(output_dir, 'final.pdb')
    
    for key in ['trajectory', 'structure']:
        path = getattr(result, key)
        if path and not os.path.exists(path):
            setattr(result, key, None)
    
    self._result = result
    return result
```

**改进建议**:
```python
def run_calvados(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
    """运行 CALVADOS 模拟
    
    建议：实际计时并在日志中输出
    """
    from .calvados_wrapper import CalvadosWrapper
    import time as time_module
    
    self._ensure_setup()
    self._ensure_not_running()
    
    dirs = self.prepare_calvados_output()
    output_dir = dirs['output_dir']
    raw_dir = dirs['raw_dir']
    
    self.is_running = True
    result = SimulationResult()
    result.output_dir = output_dir
    
    start_time = time_module.time()
    
    try:
        # 设置 GPU
        os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
        
        # 使用公共方法（如果存在）或重构为公共接口
        wrapper = CalvadosWrapper(self.config)
        # 建议：将 _write_to_dir 重命名为 write_to_dir 或提供公共接口
        files = wrapper.write(raw_dir, overwrite=True)
        
        # 运行模拟
        from multiscale2.extern.ms2_calvados.calvados import sim as calvados_sim
        calvados_sim.run(
            path=raw_dir,
            fconfig='config.yaml',
            fcomponents='components.yaml'
        )
        
        elapsed = time_module.time() - start_time
        
        # 组织输出
        wrapper._organize_output(raw_dir, output_dir, self.config.system_name)
        
        # 写入日志
        wrapper._write_simulation_log(output_dir, f"{self.config.system_name}_CG", elapsed, True)
        
        result.success = True
        print(f"  ✓ CALVADOS simulation completed ({elapsed:.1f}s)")
    
    except Exception as e:
        result.success = False
        result.errors.append(str(e))
        print(f"  ✗ CALVADOS simulation failed: {e}")
    
    finally:
        self.is_running = False
    
    # 设置结果文件路径
    result.trajectory = os.path.join(output_dir, 'trajectory.dcd')
    result.structure = os.path.join(output_dir, 'final.pdb')
    
    for key in ['trajectory', 'structure']:
        path = getattr(result, key)
        if path and not os.path.exists(path):
            setattr(result, key, None)
    
    self._result = result
    return result
```

#### 5.1.7 `_organize_calvados_output(self, raw_dir: str, output_dir: str, task_name: str)`

**功能**: 整理 CALVADOS 输出文件。

**实现逻辑**:
1. 复制轨迹文件
2. 查找并复制最终结构
3. 复制重要文件

**代码质量**: ⭐⭐⭐⭐
- 清晰的输出整理逻辑

**代码示例**:
```python
def _organize_calvados_output(self, raw_dir: str, output_dir: str, task_name: str):
    """整理 CALVADOS 输出文件到统一结构
    
    统一命名规则：
    - trajectory.dcd  <- {system_name}.dcd
    - final.pdb       <- checkpoint.pdb 或时间戳 PDB
    """
    import shutil
    
    sysname = self.config.system_name
    
    # 处理轨迹文件
    src_dcd = os.path.join(raw_dir, f'{sysname}.dcd')
    dst_dcd = os.path.join(output_dir, 'trajectory.dcd')
    if os.path.exists(src_dcd):
        shutil.copy2(src_dcd, dst_dcd)
        print(f"  📦 trajectory.dcd")
    
    # 查找最终结构
    src_pdb = os.path.join(raw_dir, 'checkpoint.pdb')
    if not os.path.exists(src_pdb):
        for f in os.listdir(raw_dir):
            if f.endswith('.pdb') and f != 'top.pdb':
                src_pdb = os.path.join(raw_dir, f)
                break
    
    dst_pdb = os.path.join(output_dir, 'final.pdb')
    if os.path.exists(src_pdb):
        shutil.copy2(src_pdb, dst_pdb)
        print(f"  📦 final.pdb")
    
    # 复制重要文件
    important_files = [
        (f'{sysname}.xml', 'system.xml'),
        ('top.pdb', 'top.pdb'),
        ('restart.chk', 'restart.chk'),
        ('checkpoint.pdb', 'checkpoint.pdb'),
    ]
    for src_name, dst_name in important_files:
        src = os.path.join(raw_dir, src_name)
        if os.path.exists(src):
            dst = os.path.join(raw_dir, dst_name)
            if src != dst:
                shutil.copy2(src, dst)
    
    print(f"  📁 原始输出已整理到: {raw_dir}")
```

#### 5.1.8 `_write_simulation_log(self, output_dir: str, task_name: str, elapsed: float, success: bool)`

**功能**: 写入模拟日志。

**代码质量**: ⭐⭐⭐⭐
- 完整的日志内容
- 格式清晰

**代码示例**:
```python
def _write_simulation_log(self, output_dir: str, task_name: str, elapsed: float, success: bool):
    """写入高层级模拟日志"""
    from datetime import datetime
    
    log_file = os.path.join(output_dir, 'simulation.log')
    
    status = "SUCCESS" if success else "FAILED"
    components_info = []
    for comp in self.config.components:
        comp_info = f"  - {comp.name}: {comp.type.value}, nmol={comp.nmol}"
        if comp.type == ComponentType.IDP:
            comp_info += f", seq={comp.ffasta}"
        elif comp.type == ComponentType.MDP:
            comp_info += f", pdb={comp.fpdb}"
        components_info.append(comp_info)
    
    log_content = f"""# Multiscale2 CG Simulation Log
# ============================

Task: {task_name}
Force Field: CALVADOS
Date: {datetime.now().isoformat()}

Status: {status}
Duration: {elapsed:.2f} seconds

System Configuration:
  Box: {self.config.box} nm
  Temperature: {self.config.temperature} K
  Ionic Strength: {self.config.ionic} M
  Topology: {self.config.topol.value}

Components ({len(self.config.components)}):
{chr(10).join(components_info)}

Output Files:
  - final.pdb: Final structure
  - trajectory.dcd: Simulation trajectory
  - raw/: Native simulation output files
"""
    with open(log_file, 'w') as f:
        f.write(log_content)
    print(f"  📝 simulation.log")
```

#### 5.1.9 `run_hps(self, gpu_id: int = 0, **kwargs) -> SimulationResult`

**功能**: 运行 HPS-Urry 模拟。

**当前状态**: 占位符实现（TODO）

**代码质量**: ⭐⭐
- 框架完整
- 缺少实际实现

**代码示例**:
```python
def run_hps(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
    """运行 HPS-Urry 模拟
    
    Args:
        gpu_id: GPU 设备 ID
        **kwargs: 额外参数
    
    Returns:
        SimulationResult
    """
    self._ensure_setup()
    self._ensure_not_running()
    
    self.is_running = True
    result = SimulationResult()
    result.output_dir = self.output_dir
    
    try:
        print(f"\n[HPS-Urry] Running simulation...")
        print(f"  GPU ID: {gpu_id}")
        
        # TODO: 实现 HPS-Urry runner
        # 使用 OpenABC 包的 HPS-Urry 力场
        
        result.success = True
        print(f"  ✓ HPS-Urry simulation completed (placeholder)")
    
    except ImportError as e:
        result.success = False
        result.errors.append(f"OpenABC not installed: {e}")
        print(f"  ✗ HPS-Urry simulation failed: OpenABC not available")
    except Exception as e:
        result.success = False
        result.errors.append(str(e))
        print(f"  ✗ HPS-Urry simulation failed: {e}")
    
    finally:
        self.is_running = False
    
    self._result = result
    return result
```

#### 5.1.10 `run_moff(self, gpu_id: int = 0, **kwargs) -> SimulationResult`

**功能**: 运行 MOFF 模拟。

**当前状态**: 占位符实现（TODO）

**代码质量**: ⭐⭐
- 框架完整
- 缺少实际实现

#### 5.1.11 `run_cocomo(self, gpu_id: int = 0, **kwargs) -> SimulationResult`

**功能**: 运行 COCOMO 模拟。

**当前状态**: 占位符实现（TODO）

**代码质量**: ⭐⭐
- 框架完整
- 缺少实际实现

#### 5.1.12 `run_openmpipi(self, gpu_id: int = 0, **kwargs) -> SimulationResult`

**功能**: 运行 OpenMpipi 模拟。

**当前状态**: 占位符实现（TODO）

**代码质量**: ⭐⭐
- 框架完整
- 缺少实际实现

#### 5.1.13 `get_result(self) -> Optional[SimulationResult]`

**功能**: 获取最近的模拟结果。

**代码质量**: ⭐⭐⭐⭐
- 简洁明了

**代码示例**:
```python
def get_result(self) -> Optional[SimulationResult]:
    """获取最近的模拟结果"""
    return self._result
```

#### 5.1.14 `cleanup(self)`

**功能**: 清理临时文件。

**当前问题**: `self.is_setup` 应该是 `self.is_setup`（属性名不匹配）

**代码质量**: ⭐⭐
- 有错误

**代码示例**:
```python
def cleanup(self):
    """清理临时文件"""
    self.is_setup = False  # 错误：应该是 self.is_setup
    self._result = None
```

**修复**:
```python
def cleanup(self):
    """清理临时文件"""
    self.is_setup = False
    self._result = None
```

---

## 六、代码审查总结

### 6.1 优点

1. **清晰的架构**: 配置、结果、模拟器职责分离
2. **完整的类型注解**: 使用 dataclass 和 Enum
3. **丰富的验证逻辑**: validate() 方法检查配置有效性
4. **序列化支持**: to_dict(), from_dict(), to_yaml(), from_yaml()
5. **统一的模拟器接口**: CGSimulator 为不同力场提供统一入口

### 6.2 问题与改进建议

#### 6.2.1 命名不一致

**问题**:
- `_DT` 和 `_FRICTION` 使用下划线前缀但不是私有属性
- `cleanup()` 中 `self.is_setup` 应该是 `self.is_setup`

**建议**: 统一命名规范

#### 6.2.2 未实现的 Runner

**问题**: `run_hps()`, `run_moff()`, `run_cocomo()`, `run_openmpipi()` 都是占位符

**建议**:
- 标记为 `@abstractmethod` 或抛出 `NotImplementedError`
- 添加明确的实现计划

#### 6.2.3 计时逻辑缺失

**问题**: `run_calvados()` 中 `elapsed` 始终为 0

**建议**: 实际计时

#### 6.2.4 类型检查不完整

**问题**: 部分方法缺少返回类型注解

**建议**: 补充类型注解

#### 6.2.5 错误处理可细化

**问题**: 错误处理过于笼统

**建议**: 添加更细粒度的错误处理

### 6.3 总体评价

| 指标 | 评分 |
|------|------|
| 代码结构 | ⭐⭐⭐⭐⭐ |
| 可读性 | ⭐⭐⭐⭐ |
| 可维护性 | ⭐⭐⭐ |
| 健壮性 | ⭐⭐⭐⭐ |
| 文档完整性 | ⭐⭐⭐⭐⭐ |

**总体评分**: ⭐⭐⭐⭐ (4/5)

---

## 七、建议的重构方案

### 7.1 统一 Runner 接口

```python
from abc import ABC, abstractmethod
from typing import Protocol

class CGRunner(Protocol):
    """CG Runner 协议"""
    
    def run(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """运行模拟"""
        ...

class BaseCGRunner(ABC):
    """CG Runner 抽象基类"""
    
    def __init__(self, config: CGSimulationConfig):
        self.config = config
        self._result: Optional[SimulationResult] = None
    
    @abstractmethod
    def run(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """运行模拟（子类必须实现）"""
        pass
    
    def get_result(self) -> Optional[SimulationResult]:
        """获取结果"""
        return self._result
```

### 7.2 计时改进

```python
def run_calvados(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
    """运行 CALVADOS 模拟（带计时）"""
    import time
    from datetime import datetime
    
    start_time = time.time()
    start_datetime = datetime.now()
    
    try:
        # ... 运行逻辑 ...
        elapsed = time.time() - start_time
        
    except Exception as e:
        elapsed = time.time() - start_time
        # ... 错误处理 ...
    
    # 在日志中使用 start_datetime 记录开始时间
    return result
```

### 7.3 错误处理细化

```python
class CGSimulatorError(Exception):
    """CG 模拟器基础异常"""
    pass

class ConfigurationError(CGSimulatorError):
    """配置错误"""
    pass

class SimulationSetupError(CGSimulatorError):
    """设置错误"""
    pass

class SimulationRuntimeError(CGSimulatorError):
    """运行时错误"""
    pass
```

---

## 八、测试建议

### 8.1 单元测试

```python
import pytest
from multiscale2.src.cg import (
    CGSimulationConfig,
    CGComponent,
    CGSimulator,
    ComponentType,
    TopologyType,
    Platform,
    SimulationParams,
    SimulationResult,
)

class TestCGSimulationConfig:
    """CGSimulationConfig 测试类"""
    
    def test_default_values(self):
        """测试默认值"""
        config = CGSimulationConfig()
        
        assert config.system_name == "cg_simulation"
        assert config.box == [25.0, 25.0, 30.0]
        assert config.temperature == 310.0
        assert config.ionic == 0.15
        assert config.topol == TopologyType.CUBIC
        assert config.output_dir == "output_cg"
    
    def test_add_component(self):
        """测试添加组件"""
        config = CGSimulationConfig()
        comp = CGComponent(
            name="test",
            type=ComponentType.IDP,
            nmol=5,
        )
        
        config.add_component(comp)
        
        assert len(config.components) == 1
        assert config.total_molecules() == 5
    
    def test_validate_empty_config(self):
        """测试空配置验证"""
        config = CGSimulationConfig()
        errors = config.validate()
        
        assert "At least one component is required" in errors
    
    def test_validate_missing_fasta(self):
        """测试缺少 FASTA 验证"""
        config = CGSimulationConfig()
        config.add_component(CGComponent(
            name="test",
            type=ComponentType.IDP,
            nmol=1,
            # ffasta 未设置
        ))
        
        errors = config.validate()
        
        assert any("IDP requires ffasta file" in e for e in errors)
    
    def test_validate_missing_pdb(self):
        """测试缺少 PDB 验证"""
        config = CGSimulationConfig()
        config.add_component(CGComponent(
            name="test",
            type=ComponentType.MDP,
            nmol=1,
            restraint=True,
            # fpdb 未设置
        ))
        
        errors = config.validate()
        
        assert any("MDP requires fpdb file" in e for e in errors)
    
    def test_to_dict_and_from_dict(self):
        """测试序列化和反序列化"""
        config = CGSimulationConfig(
            system_name="test_sys",
            box=[30.0, 30.0, 40.0],
            temperature=298.0,
        )
        
        config.add_component(CGComponent(
            name="prot1",
            type=ComponentType.IDP,
            nmol=10,
        ))
        
        # 序列化
        d = config.to_dict()
        
        # 反序列化
        config2 = CGSimulationConfig.from_dict(d)
        
        assert config2.system_name == "test_sys"
        assert config2.box == [30.0, 30.0, 40.0]
        assert config2.temperature == 298.0
        assert len(config2.components) == 1


class TestCGSimulator:
    """CGSimulator 测试类"""
    
    @pytest.fixture
    def sample_config(self):
        """示例配置"""
        config = CGSimulationConfig(
            system_name="test_sim",
        )
        config.add_component(CGComponent(
            name="test_protein",
            type=ComponentType.IDP,
            nmol=5,
            ffasta="tests/test_data/test.fasta",
        ))
        return config
    
    def test_init_valid_config(self, sample_config):
        """测试有效配置初始化"""
        sim = CGSimulator(sample_config)
        
        assert sim.is_setup is False
        assert sim.is_running is False
    
    def test_init_invalid_config(self):
        """测试无效配置初始化"""
        config = CGSimulationConfig()
        
        with pytest.raises(ValueError):
            CGSimulator(config)
    
    def test_setup(self, sample_config, tmp_path):
        """测试设置"""
        sim = CGSimulator(sample_config)
        
        output_dir = str(tmp_path / "output")
        result = sim.setup(output_dir)
        
        assert sim.is_setup is True
        assert os.path.exists(result['output_dir'])
        assert os.path.exists(result['config'])
    
    def test_setup_already_exists(self, sample_config, tmp_path):
        """测试设置目录已存在"""
        sim = CGSimulator(sample_config)
        
        output_dir = str(tmp_path / "output")
        sim.setup(output_dir)
        
        with pytest.raises(FileExistsError):
            sim.setup(output_dir, overwrite=False)
    
    def test_ensure_setup(self, sample_config):
        """测试确保设置"""
        sim = CGSimulator(sample_config)
        
        with pytest.raises(RuntimeError):
            sim._ensure_setup()
    
    def test_ensure_not_running(self, sample_config):
        """测试确保未运行"""
        sim = CGSimulator(sample_config)
        
        sim.is_running = True
        
        with pytest.raises(RuntimeError):
            sim._ensure_not_running()
    
    def test_cleanup(self, sample_config):
        """测试清理"""
        sim = CGSimulator(sample_config)
        
        sim.is_setup = True
        sim._result = SimulationResult()
        
        sim.cleanup()
        
        assert sim.is_setup is False
        assert sim._result is None
```

---

## 九、API 参考

### 9.1 CGSimulationConfig

```python
class CGSimulationConfig:
    system_name: str                    # 系统名称
    box: List[float]                    # 盒子尺寸 [x, y, z]
    temperature: float                  # 温度 (K)
    ionic: float                        # 离子强度 (M)
    topol: TopologyType                 # 拓扑类型
    simulation: SimulationParams        # 模拟参数
    components: List[CGComponent]       # 组件列表
    output_dir: str                     # 输出目录
    
    def add_component(component: CGComponent)
    def get_component(name: str) -> Optional[CGComponent]
    def total_molecules() -> int
    def validate() -> List[str]
    def to_dict() -> Dict
    def to_yaml(path: str = None)
    @classmethod from_dict(d: Dict) -> CGSimulationConfig
    @classmethod from_yaml(path: str) -> CGSimulationConfig
```

### 9.2 CGComponent

```python
class CGComponent:
    name: str                           # 组件名称
    type: ComponentType                 # 组件类型
    nmol: int                           # 分子数
    ffasta: Optional[str]               # FASTA 文件
    fpdb: Optional[str]                 # PDB 文件
    fdomains: Optional[str]             # 域定义文件
    fpae: Optional[str]                 # PAE 文件
    restraint: bool                     # 是否约束
    restraint_type: str                 # 约束类型
    use_com: bool                       # 使用质心
    k_harmonic: float                   # 谐波力常数
    colabfold: int                      # Colabfold 格式
    charge_termini: str                 # 末端电荷
    
    def validate() -> List[str]
    def to_dict() -> Dict
    @classmethod from_dict(d: Dict) -> CGComponent
```

### 9.3 CGSimulator

```python
class CGSimulator:
    config: CGSimulationConfig
    output_dir: Optional[str]
    is_setup: bool
    is_running: bool
    _result: Optional[SimulationResult]
    
    def __init__(config: CGSimulationConfig)
    def setup(output_dir: str, overwrite: bool = False) -> Dict[str, str]
    def prepare_calvados_output() -> Dict[str, str]
    def run_calvados(gpu_id: int = 0, **kwargs) -> SimulationResult
    def run_hps(gpu_id: int = 0, **kwargs) -> SimulationResult
    def run_moff(gpu_id: int = 0, **kwargs) -> SimulationResult
    def run_cocomo(gpu_id: int = 0, **kwargs) -> SimulationResult
    def run_openmpipi(gpu_id: int = 0, **kwargs) -> SimulationResult
    def get_result() -> Optional[SimulationResult]
    def cleanup()
```

---

**报告完成**

*审查人: AI Assistant*  
*审查日期: 2024年12月30日*

