#!/usr/bin/env python3
"""
Coarse-Grained Simulation Module

A unified interface for coarse-grained molecular dynamics simulations
supporting multiple force fields (CALVADOS, HPS, MOFF, COCOMO, OpenMpipi).

Architecture:
- CGSimulationConfig: Configuration data class
- CGComponent: Individual component specification
- CGSimulator: Main simulator class with runner methods

Usage:
    from multiscale2.src import CGSimulationConfig, CGSimulator
    
    config = CGSimulationConfig.from_yaml("config.yaml")
    sim = CGSimulator(config)
    sim.setup("output/")
    sim.run_calvados(gpu_id=0)
"""

import os
import yaml
import shutil
import warnings
import numpy as np
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Any
from enum import Enum
from pathlib import Path

from openmm import Platform, LangevinIntegrator, LangevinMiddleIntegrator, XmlSerializer
import openmm as mm
import openmm.unit as unit
from openmm.app import (
    Simulation, StateDataReporter, PDBFile
)
from openmm.unit import (
    picoseconds, nanometers, kilojoule, mole,
    Quantity, kilocalorie, amu, kelvin
)

# mdtraj reporters for trajectory saving with PBC support
from mdtraj.reporters import XTCReporter

from .pdb_tool import ChainLabel


# =============================================================================
# Enums
# =============================================================================

class ComponentType(Enum):
    """组件类型"""
    IDP = "idp"   # 无序蛋白 - 基于序列
    MDP = "mdp"   # 折叠蛋白 - 基于结构


class TopologyType(Enum):
    """拓扑类型"""
    CUBIC = "cubic"   # 立方体盒子（网格放置）
    SLAB = "slab"     # 平面限制（相分离）


class ComputePlatform(Enum):
    """计算平台"""
    CPU = "CPU"
    CUDA = "CUDA"


# =============================================================================
# Configuration Classes
# =============================================================================

@dataclass
class SimulationParams:
    """
    核心模拟参数

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
    platform: ComputePlatform = ComputePlatform.CUDA
    verbose: bool = True

    def to_dict(self) -> Dict:
        d = {
            'steps': self.steps,
            'wfreq': self.wfreq,
            'platform': self.platform.value,
            'verbose': self.verbose,
        }
        # 移除 None 值
        return {k: v for k, v in d.items() if v is not None}

    @classmethod
    def from_dict(cls, d: Dict) -> 'SimulationParams':
        if 'platform' in d and isinstance(d['platform'], str):
            d['platform'] = ComputePlatform(d['platform'])
        # 移除 None 值
        d = {k: v for k, v in d.items() if v is not None}
        return cls(**d)


@dataclass
class CGComponent:
    """
    单个组件规格
    
    Attributes:
        name: 组件唯一标识
        type: 组件类型（IDP 或 MDP）
        nmol: 该类型的分子数
        ffasta: FASTA 文件路径（IDP）
        fpdb: PDB 文件路径（MDP）
        fdomains: 域定义文件路径（MDP）
        restraint: 是否应用结构约束（MDP）
        restraint_type: 约束类型（harmonic, go）
        use_com: 是否在质心应用约束
        k_harmonic: 谐波约束力常数
        colabfold: PAE 格式（0=EBI, 1&2=Colabfold）
        charge_termini: 末端电荷（both, n, c, none）
    """
    name: str
    type: ComponentType
    nmol: int = 1
    
    # 输入文件
    ffasta: Optional[str] = None       # FASTA 文件（IDP）
    fpdb: Optional[str] = None         # PDB 文件（MDP）
    fdomains: Optional[str] = None     # 域定义文件（MDP）
    fpae: Optional[str] = None         # PAE JSON（Go 势）
    
    # 约束设置
    restraint: bool = False
    restraint_type: str = "harmonic"
    use_com: bool = False
    k_harmonic: float = 700.0
    colabfold: int = 1
    
    # 电荷设置
    charge_termini: str = "both"
    
    # 派生属性
    seq: Optional[str] = None
    nres: int = 0
    
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
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'CGComponent':
        comp_type = d.get('type', 'idp')
        if isinstance(comp_type, str):
            # 转换为小写以匹配枚举值
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
    
    def validate(self) -> List[str]:
        """验证配置
        
        fdomains 支持两种格式：
        1. 文件路径：'domains.yaml' - 检查文件是否存在
        2. 内联 YAML：'TDP43:\n  - [3, 76]\n...' - 不检查
        """
        errors = []
        
        def _is_inline_yaml(text: str) -> bool:
            """检查是否是内联 YAML（而不是文件路径）"""
            if not text:
                return False
            stripped = text.strip()
            # 以 {、[、字母开头，或者包含换行符且有 YAML 特征
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
                # 内联 YAML 不需要检查文件存在
                if not _is_inline_yaml(self.fdomains) and not os.path.exists(self.fdomains):
                    errors.append(f"Component '{self.name}': Domains file not found: {self.fdomains}")
        return errors


@dataclass
class CGSimulationConfig:
    """
    完整模拟配置
    
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
    ionic: float = 0.15              # Molar (离子强度)
    
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
    
    def add_component(self, component: CGComponent):
        """添加组件"""
        self.components.append(component)
    
    def get_component(self, name: str) -> Optional[CGComponent]:
        """根据名称获取组件"""
        for comp in self.components:
            if comp.name == name:
                return comp
        return None
    
    def total_molecules(self) -> int:
        """计算总分子数"""
        return sum(comp.nmol for comp in self.components)
    
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
    
    def to_yaml(self, path: str = None):
        """保存到 YAML 文件"""
        d = self.to_dict()
        if path:
            with open(path, 'w') as f:
                yaml.dump(d, f, default_flow_style=False, sort_keys=False)
        else:
            return yaml.dump(d, default_flow_style=False, sort_keys=False)
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'CGSimulationConfig':
        """从字典创建"""
        topol = d.get('topol', 'cubic')
        if isinstance(topol, str):
            # 转换为小写以匹配枚举值
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
    
    @classmethod
    def from_yaml(cls, path: str) -> 'CGSimulationConfig':
        """从 YAML 文件加载

        自动将相对路径转换为相对于配置文件的绝对路径
        """
        import os

        # 获取配置文件的目录，用于解析相对路径
        config_dir = os.path.dirname(os.path.abspath(path))

        with open(path, 'r', encoding='utf-8') as f:
            d = yaml.safe_load(f)
        d['config_path'] = path

        config = cls.from_dict(d)

        # 转换组件中的相对路径为绝对路径
        for comp in config.components:
            # 处理 ffasta 路径
            if comp.ffasta and not os.path.isabs(comp.ffasta):
                comp.ffasta = os.path.join(config_dir, comp.ffasta)

            # 处理 fpdb 路径
            if comp.fpdb and not os.path.isabs(comp.fpdb):
                comp.fpdb = os.path.join(config_dir, comp.fpdb)

            # 处理 fdomains 路径（跳过内联 YAML）
            if comp.fdomains and not cls._is_inline_yaml(comp.fdomains):
                if not os.path.isabs(comp.fdomains):
                    comp.fdomains = os.path.join(config_dir, comp.fdomains)

            # 处理 fpae 路径
            if comp.fpae and not os.path.isabs(comp.fpae):
                comp.fpae = os.path.join(config_dir, comp.fpae)

        return config

    @staticmethod
    def _is_inline_yaml(text: str) -> bool:
        """检查是否是内联 YAML（而不是文件路径）"""
        if not text:
            return False
        stripped = text.strip()
        # 以 {、[开头，或者包含换行符且有 YAML 特征
        if stripped.startswith('{') or stripped.startswith('['):
            return True
        if '\n' in stripped and (':' in stripped or stripped.startswith('-')):
            return True
        return False


# =============================================================================
# Simulation Result
# =============================================================================

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


@dataclass
class TopologyInfo:
    """
    拓扑信息缓存

    存储体系的完整拓扑信息，用于各种分析接口。

    Attributes:
        global_sequence: 全局序列（所有残基拼接）
        chain_ids: 每个残基的链ID（纯数字，从1开始）
        is_folded: 每个残基是否属于folded domain（0=unfolded/IDP，1=folded domain）
        molecule_indices: 每个残基属于哪个分子 (1, 2, 3, ...)
        component_names: 每个残基属于哪个组件
        local_residue_indices: 在单链中的残基编号 (1-based)
        sasa_values: 每个残基的SASA值（从全原子结构计算）
    """
    global_sequence: str = ""
    chain_ids: List[int] = field(default_factory=list)
    is_folded: List[int] = field(default_factory=list)
    molecule_indices: List[int] = field(default_factory=list)
    component_names: List[str] = field(default_factory=list)
    local_residue_indices: List[int] = field(default_factory=list)
    sasa_values: List[float] = field(default_factory=list)


# =============================================================================
# CG Simulator with Multiple Runners
# =============================================================================

class CGSimulator:
    """
    粗粒化模拟器
    
    支持多种力场的统一接口，包含多个 runner 方法。
    
    Attributes:
        config: 模拟配置
        output_dir: 输出目录
        is_setup: 是否已完成设置
        is_running: 是否正在运行
    """
    
    def __init__(self, config: CGSimulationConfig):
        """
        初始化模拟器

        Args:
            config: CGSimulationConfig 实例
        """
        self.config = config
        self.output_dir: Optional[str] = None
        self.is_setup: bool = False
        self.is_running: bool = False
        self._result: Optional[SimulationResult] = None
        self._topology_info: Optional[TopologyInfo] = None  # 拓扑信息缓存

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
    
    # -------------------------------------------------------------------------
    # Setup Methods
    # -------------------------------------------------------------------------

    def setup(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
        """
        设置模拟环境（通用准备）

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
            # 如果 overwrite=True，不删除目录，让后续的力场特定方法来处理备份逻辑

        self.output_dir = output_dir

        print(f"\n[CGSimulator] Setting up...")
        print(f"  Output directory: {output_dir}")
        print(f"  System: {self.config.system_name}")

        self.is_setup = True
        print(f"  ✓ Setup complete")

        return {
            'output_dir': output_dir,
            'config': os.path.join(output_dir, 'config.yaml'),
        }

    # -------------------------------------------------------------------------
    # Topology / Sequence Interface Methods
    # -------------------------------------------------------------------------

    def get_composition(self) -> List[Dict[str, Any]]:
        """
        返回整个体系的构成信息

        Returns:
            List of dicts, each containing:
            - name: 组件名称
            - nmol: 分子数
            - sequence: 单链序列
            - nres: 序列长度
            - type: 组件类型 (IDP 或 MDP)
        """
        composition = []
        for comp in self.config.components:
            # 获取序列
            seq = self._get_component_sequence(comp)
            composition.append({
                'name': comp.name,
                'nmol': comp.nmol,
                'sequence': seq,
                'nres': len(seq),
                'type': comp.type.value,
            })
        return composition

    def get_global_sequence(self) -> str:
        """
        返回全局序列（将所有组分的所有分子按顺序拼接）

        Returns:
            全局序列字符串
        """
        info = self._get_topology_info()
        return info.global_sequence

    def get_chain_ids(self) -> List[int]:
        """
        返回与全局序列等长的链ID列表

        链ID为纯数字，从1开始编号。

        Returns:
            链ID列表（纯数字，从1开始）
        """
        info = self._get_topology_info()
        return info.chain_ids

    def get_folded_domains(self) -> List[int]:
        """
        返回与全局序列等长的folded domain标记列表

        Returns:
            整数列表，0 表示 IDP 或 unfolded 区域，1 表示 folded domain
        """
        info = self._get_topology_info()
        return info.is_folded

    def get_chain_identifiers(self) -> List[str]:
        """
        返回与全局序列等长的链标识符列表

        返回适合用于 OpenMM topology 的链标识符字符串。
        格式：'A', 'B', ..., 'Z', 'A1', 'B1', ..., 'Z1', 'A2', ...
        支持无限多的链。

        Returns:
            链标识符列表，长度与全局序列相同
        """
        info = self._get_topology_info()
        chain_ids = info.chain_ids

        def chain_id_to_identifier(chain_id: int) -> str:
            """Convert chain ID to unique identifier."""
            if chain_id <= 26:
                return chr(ord('A') + chain_id - 1)
            else:
                letter_idx = (chain_id - 1) % 26
                suffix = (chain_id - 1) // 26
                return chr(ord('A') + letter_idx) + str(suffix)

        return [chain_id_to_identifier(cid) for cid in chain_ids]

    def get_unique_chain_identifiers(self) -> List[str]:
        """
        返回所有唯一的链标识符列表

        Returns:
            唯一链标识符列表，按字母顺序排序
        """
        identifiers = set(self.get_chain_identifiers())
        return sorted(list(identifiers))

    def _build_component_names(self) -> List[str]:
        """
        构建每个残基对应的组件名称列表

        Returns:
            组件名称列表，长度与全局序列相同
        """
        comp_names = []
        for comp in self.config.components:
            for _ in range(comp.nmol):
                comp_names.extend([comp.name] * len(self._get_component_sequence(comp)))
        return comp_names

    def _compute_sasa_for_component(self, comp: CGComponent) -> List[float]:
        """
        为 MDP 组件计算 SASA 值

        使用 mdsim 库从全原子 PDB 文件计算每个残基的 SASA。
        对于 IDP 组件，返回默认 SASA 值。

        Args:
            comp: CGComponent 实例

        Returns:
            SASA 值列表（每个残基一个值）
        """
        # IDP 使用默认 SASA 值
        if comp.type == ComponentType.IDP:
            nres = len(self._get_component_sequence(comp))
            return [5.0] * nres

        # MDP 需要全原子 PDB 文件
        if comp.type == ComponentType.MDP:
            if not comp.fpdb:
                raise ValueError(f"Component '{comp.name}' is MDP but no fpdb file specified")

            # 尝试导入 system_handling 模块
            try:
                from multiscale2.extern.cocomo.src.cocomo.system_handling import ComponentType as SASAComponentType
            except ImportError:
                warnings.warn(
                    f"Could not import cocomo.system_handling for SASA calculation. "
                    f"Using default SASA values for component '{comp.name}'.",
                    UserWarning
                )
                nres = len(self._get_component_sequence(comp))
                return [5.0] * nres

            try:
                # 创建 ComponentType 并自动计算 SASA
                ctype = SASAComponentType(
                    name=comp.name,
                    pdb=comp.fpdb,
                    getsasa="auto",  # 使用 mdsim 自动计算
                    mask_sasa_bydomain=True  # 按 domain 掩码
                )

                # 提取 SASA 值列表
                sasa_list = [val for _, val in ctype.sasa]

                # 检查是否成功计算了 SASA
                nres = len(self._get_component_sequence(comp))
                if len(sasa_list) != nres:
                    warnings.warn(
                        f"SASA length ({len(sasa_list)}) doesn't match sequence length ({nres}) "
                        f"for component '{comp.name}'. Padding with default values.",
                        UserWarning
                    )
                    # 填充或截断到正确的长度
                    if len(sasa_list) < nres:
                        sasa_list.extend([5.0] * (nres - len(sasa_list)))
                    else:
                        sasa_list = sasa_list[:nres]

                return sasa_list

            except Exception as e:
                warnings.warn(
                    f"Failed to compute SASA for component '{comp.name}': {e}. "
                    f"Using default SASA values.",
                    UserWarning
                )
                nres = len(self._get_component_sequence(comp))
                return [5.0] * nres

        # 默认情况
        nres = len(self._get_component_sequence(comp))
        return [5.0] * nres

    def _compute_all_sasa_values(self) -> List[float]:
        """
        计算整个系统的所有 SASA 值

        为每个组件计算 SASA 值，并根据分子数量重复。

        Returns:
            SASA 值列表，长度与全局序列相同
        """
        all_sasa = []

        for comp in self.config.components:
            # 计算单链的 SASA
            single_sasa = self._compute_sasa_for_component(comp)

            # 根据分子数重复
            for _ in range(comp.nmol):
                all_sasa.extend(single_sasa)

        return all_sasa

    def _get_sasa_values(self) -> Optional[np.ndarray]:
        """
        尝试加载或计算 SASA 值

        优先从文件加载，如果文件不存在则计算。
        对于 MDP 组件，使用 mdsim 从全原子 PDB 计算。
        对于 IDP 组件，使用默认 SASA 值。

        Returns:
            SASA 值数组，如果无法计算则返回 None
        """
        # 尝试加载缓存文件
        surface_file = os.path.join(self.output_dir or self.config.output_dir, 'surface')
        if os.path.exists(surface_file):
            print(f"  加载 SASA 数据: {surface_file}")
            return np.loadtxt(surface_file)

        # 尝试其他可能的路径
        alt_surface_file = os.path.join(self.output_dir or self.config.output_dir, 'sasa_values.txt')
        if os.path.exists(alt_surface_file):
            print(f"  加载 SASA 数据: {alt_surface_file}")
            return np.loadtxt(alt_surface_file)

        # 如果没有缓存文件，尝试计算
        try:
            sasa = self._compute_all_sasa_values()
            if sasa:
                print(f"  计算得到 {len(sasa)} 个 SASA 值")
                return np.array(sasa)
        except Exception as e:
            print(f"  SASA 计算失败: {e}")

        print(f"  未找到 SASA 文件，使用默认值")
        return None

    def get_sasa_values(self) -> List[float]:
        """
        获取整个系统的 SASA 值列表

        返回与全局序列等长的 SASA 值列表。
        对于 MDP 组件，从全原子 PDB 文件计算；
        对于 IDP 组件，使用默认 SASA 值。

        Returns:
            SASA 值列表
        """
        info = self._get_topology_info()
        if info.sasa_values:
            return info.sasa_values

        # 如果拓扑信息中没有，计算
        sasa = self._compute_all_sasa_values()

        # 缓存到拓扑信息
        self._topology_info.sasa_values = sasa

        return sasa

    def _get_topology_info(self) -> TopologyInfo:
        """
        获取拓扑信息（懒加载，缓存计算结果）

        Returns:
            TopologyInfo 实例
        """
        if self._topology_info is None:
            self._topology_info = self._build_topology_info()
        return self._topology_info

    def _build_topology_info(self) -> TopologyInfo:
        """
        构建拓扑信息

        构建全局序列、链ID、folded domain等信息。
        对于 MDP 组件，从全原子 PDB 文件计算 SASA 值。

        Returns:
            TopologyInfo 实例
        """
        global_sequence_parts = []
        chain_ids = []
        is_folded = []
        molecule_indices = []
        component_names = []
        local_residue_indices = []
        sasa_values = []

        # 当前链ID（从1开始）
        current_chain_id = 1

        # 预计算每个组件的 SASA 值（只计算一次）
        component_sasa = {}
        for comp in self.config.components:
            component_sasa[comp.name] = self._compute_sasa_for_component(comp)

        for comp in self.config.components:
            # 获取单链序列
            single_seq = self._get_component_sequence(comp)
            nres = len(single_seq)
            nmol = comp.nmol

            # 获取单链的folded domain信息
            single_folded = self._get_component_folded_domains(comp, nres)

            # 获取该组件的 SASA 值
            single_sasa = component_sasa.get(comp.name, [5.0] * nres)

            # 为每个分子构建信息
            for mol_idx in range(nmol):
                # 添加序列
                global_sequence_parts.append(single_seq)

                # 添加每个残基的信息
                for res_idx in range(nres):
                    chain_ids.append(current_chain_id)
                    is_folded.append(single_folded[res_idx])
                    molecule_indices.append(current_chain_id)  # 链ID就是分子ID
                    component_names.append(comp.name)
                    local_residue_indices.append(res_idx + 1)  # 1-based
                    sasa_values.append(single_sasa[res_idx])

                current_chain_id += 1

        return TopologyInfo(
            global_sequence="".join(global_sequence_parts),
            chain_ids=chain_ids,
            is_folded=is_folded,
            molecule_indices=molecule_indices,
            component_names=component_names,
            local_residue_indices=local_residue_indices,
            sasa_values=sasa_values,
        )

    def _get_component_sequence(self, comp: CGComponent) -> str:
        """
        获取组件的单链序列

        Args:
            comp: CGComponent 实例

        Returns:
            单链序列字符串
        """
        # 如果已有序列，直接返回
        if comp.seq:
            return comp.seq

        # 从FASTA文件读取（IDP）或PDB文件读取（MDP）
        if comp.type == ComponentType.IDP:
            if comp.ffasta:
                return self._read_fasta(comp.ffasta, component_name=comp.name)
            else:
                raise ValueError(f"Component '{comp.name}' is IDP but no ffasta file specified")
        elif comp.type == ComponentType.MDP:
            if comp.fpdb:
                return self._seq_from_pdb(comp.fpdb)
            else:
                raise ValueError(f"Component '{comp.name}' is MDP but no fpdb file specified")

        return ""

    def _read_fasta(self, fasta_path: str, component_name: str = None) -> str:
        """
        从FASTA文件读取序列

        如果指定了 component_name，则返回匹配该名称的序列。
        如果未指定或找不到匹配，则返回第一条序列。

        Args:
            fasta_path: FASTA文件路径
            component_name: 组件名称（用于选择序列）

        Returns:
            序列字符串
        """
        from Bio import SeqIO

        # 处理相对路径
        if not os.path.isabs(fasta_path):
            # 相对于当前工作目录或配置文件目录
            if hasattr(self.config, 'config_path') and self.config.config_path:
                config_dir = os.path.dirname(os.path.abspath(self.config.config_path))
                fasta_path = os.path.join(config_dir, fasta_path)

        # 读取FASTA文件
        records = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            raise ValueError(f"Empty or invalid FASTA file: {fasta_path}")

        # 如果指定了组件名称，尝试匹配
        if component_name:
            if component_name in records:
                return str(records[component_name].seq)
            else:
                # 尝试不区分大小写的匹配
                for name in records:
                    if name.lower() == component_name.lower():
                        return str(records[name].seq)
                # 找不到匹配，使用第一条序列并警告
                print(f"  [WARNING] Component '{component_name}' not found in fasta, using first sequence")

        # 返回第一条序列
        return str(list(records.values())[0].seq)

    def _seq_from_pdb(self, pdb_path: str) -> str:
        """
        从PDB文件提取序列

        Args:
            pdb_path: PDB文件路径

        Returns:
            序列字符串
        """
        # 处理相对路径
        if not os.path.isabs(pdb_path):
            # 相对于当前工作目录或配置文件目录
            if hasattr(self.config, 'config_path') and self.config.config_path:
                config_dir = os.path.dirname(os.path.abspath(self.config.config_path))
                pdb_path = os.path.join(config_dir, pdb_path)

        # 使用MDAnalysis提取序列
        try:
            from MDAnalysis import Universe
        except ImportError:
            raise ImportError("MDAnalysis is required to extract sequences from PDB files")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = Universe(pdb_path)

            # 获取唯一残基
            residues = u.residues
            n_res = len(residues)

            # 3-letter to 1-letter amino acid mapping
            aa_3to1 = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                'SEC': 'U', 'PYL': 'O',  # Selenocysteine, Pyrrolysine
            }

            fastapdb = ""
            for res in residues:
                resname = res.resname
                fastapdb += aa_3to1.get(resname, 'X')

            return fastapdb

    def _get_component_folded_domains(self, comp: CGComponent, nres: int) -> List[int]:
        """
        获取组件的folded domain信息

        Args:
            comp: CGComponent 实例
            nres: 序列长度

        Returns:
            长度为nres的列表，0=unfolded/IDP，1=folded domain
        """
        # IDP默认全是0
        if comp.type == ComponentType.IDP:
            return [0] * nres

        # MDP检查是否有fdomains配置
        if not comp.fdomains:
            return [0] * nres

        # 解析fdomains YAML
        domains = self._parse_fdomains(comp.fdomains)

        # 构建folded数组
        folded = [0] * nres
        for (start, end) in domains:
            # 确保在有效范围内
            start = max(1, start)  # 1-based
            end = min(nres, end)
            if start <= end:
                for i in range(start - 1, end):  # 转换为0-based
                    folded[i] = 1

        return folded

    def _parse_fdomains(self, fdomains: str) -> List[tuple]:
        """
        解析fdomains配置

        支持两种格式：
        1. 文件路径：解析YAML文件
        2. 内联YAML：直接解析字符串

        Args:
            fdomains: fdomains配置

        Returns:
            域列表，每个域为(start, end)元组（1-based）
        """
        # 检查是否是内联YAML
        if self.config._is_inline_yaml(fdomains):
            # 直接解析字符串
            data = yaml.safe_load(fdomains)
        else:
            # 解析文件
            fdomains_abs = fdomains
            if not os.path.isabs(fdomains):
                if hasattr(self.config, 'config_path') and self.config.config_path:
                    config_dir = os.path.dirname(os.path.abspath(self.config.config_path))
                    fdomains_abs = os.path.join(config_dir, fdomains)

            if not os.path.exists(fdomains_abs):
                return []

            with open(fdomains_abs, 'r') as f:
                data = yaml.safe_load(f)

        # 解析域定义
        domains = []
        if isinstance(data, dict):
            for protein_name, domain_list in data.items():
                if isinstance(domain_list, list):
                    for domain in domain_list:
                        if isinstance(domain, (list, tuple)) and len(domain) == 2:
                            domains.append((domain[0], domain[1]))

        return domains

    def clear_topology_cache(self):
        """
        清除拓扑信息缓存

        下次调用接口方法时会重新计算。
        """
        self._topology_info = None

    def prepare_calvados_output(self) -> Dict[str, str]:
        """
        准备 CALVADOS 输出的目录结构

        统一输出结构：
        {output_dir}/
        ├── {system_name}_CG/
        │   ├── raw/                  # 原生输出
        │   ├── trajectory.xtc        # 整理后的轨迹
        │   ├── final.pdb             # 整理后的最终结构
        │   └── simulation.log        # 高层级日志

        Returns:
            包含输出路径的字典
        """
        self._ensure_setup()
        self._ensure_not_running()

        # 检查 self.output_dir 是否已包含 _CG 后缀
        expected_suffix = f"{self.config.system_name}_CG"
        if self.output_dir.endswith(expected_suffix):
            # 已包含 _CG 后缀，直接使用
            output_dir = self.output_dir
            task_name = expected_suffix
        else:
            # 添加 _CG 后缀
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

        # 先创建父目录，再创建 raw 子目录
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(raw_dir, exist_ok=True)

        return {
            'output_dir': output_dir,
            'raw_dir': raw_dir,
            'task_name': task_name,
        }
    
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
    
    def _ensure_setup(self):
        """确保已完成设置"""
        if not self.is_setup:
            raise RuntimeError("Simulation not set up. Call setup() first.")
    
    def _ensure_not_running(self):
        """确保未在运行"""
        if self.is_running:
            raise RuntimeError("Simulation is already running")

    # -------------------------------------------------------------------------
    # Pre-equilibration Methods (使用 CALVADOS 构建初始结构)
    # -------------------------------------------------------------------------
    # 说明：所有非 CALVADOS 的 runner 都会自动调用此方法进行预平衡
    # 预平衡参数是硬编码的，每个 runner 可以有不同的默认参数
    # -------------------------------------------------------------------------

    def _run_pre_equilibration(
        self,
        gpu_id: int = 0,
        steps: int = 100000,
        mapping: str = "ca",
        k_restraint: float = 10000.0,
        use_com: bool = True,
        platform: ComputePlatform = ComputePlatform.CUDA,
    ) -> Optional[str]:
        """
        运行预平衡（使用 CALVADOS 构建初始结构）

        预平衡的目的：
        1. 用 CALVADOS 力场构建初始 CG 结构
        2. 对 MDP 蛋白进行短时间模拟（带约束），使其适应 CG 表示
        3. 对 IDP 蛋白进行短时间模拟（无约束），生成初始 CG 结构
        4. 为后续的目标力场（COCOMO/HPS/MOFF/OpenMpipi）提供良好的初始结构

        对于 MDP 蛋白的映射选择：
        - CA (alpha carbon): 使用 alpha-carbon 坐标，适用于需要保留骨架信息的场景
        - COM (center of mass): 使用残基质心，适用于更平滑的映射

        Args:
            gpu_id: GPU 设备 ID
            steps: 预平衡步数（默认 10w 步）
            mapping: 映射方式 ('ca' 或 'com')
            k_restraint: 约束力常数 (kJ/(mol·nm²))
            use_com: 是否使用 COM 约束
            platform: 计算平台（CUDA 或 CPU）

        Returns:
            预平衡后的 final.pdb 文件路径，如果没有组件则返回 None
        """
        from .calvados_wrapper import CalvadosWrapper
        import shutil

        self._ensure_setup()
        self._ensure_not_running()

        # 检查是否有组件需要处理
        has_components = len(self.config.components) > 0
        has_mdp = any(comp.type == ComponentType.MDP for comp in self.config.components)
        if not has_components:
            return None

        # 对于 MDP 系统打印约束信息，IDP 系统不打印
        if has_mdp:
            print(f"\n[Pre-equilibration] 使用 CALVADOS 构建初始结构（MDP 蛋白，带约束）...")
            print(f"  映射方式: {mapping.upper()}")
            print(f"  约束力常数: {k_restraint} kJ/(mol·nm²)")
            print(f"  使用 COM 约束: {'是' if use_com else '否'}")
        else:
            print(f"\n[Pre-equilibration] 使用 CALVADOS 构建初始结构（纯 IDP 系统）...")
            print(f"  映射方式: {mapping.upper()}")
        print(f"  预平衡步数: {steps}")

        # 保存原始组件配置
        original_components = []
        for comp in self.config.components:
            original_components.append({
                'restraint': comp.restraint,
                'restraint_type': comp.restraint_type,
                'use_com': comp.use_com,
                'k_harmonic': comp.k_harmonic,
            })

        # 创建临时配置用于预平衡
        temp_config = self._create_preequil_config(
            steps=steps,
            mapping=mapping,
            k_restraint=k_restraint,
            use_com=use_com,
            platform=platform,
        )

        # 输出目录结构：{output_dir}/{system_name}_CG/equilibration/
        # 与 run_calvados 保持一致
        task_name = f"{self.config.system_name}_CG"
        equilibration_dir = os.path.join(self.output_dir, task_name, 'equilibration')
        raw_dir = os.path.join(equilibration_dir, 'raw')

        # 如果目录已存在，先备份
        if os.path.exists(equilibration_dir):
            import time
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{equilibration_dir}_backup_{timestamp}"
            shutil.move(equilibration_dir, backup_dir)
            print(f"  📁 备份旧预平衡结果到: {backup_dir}")

        os.makedirs(raw_dir, exist_ok=True)

        try:
            # 根据平台设置 GPU 环境变量
            if platform == ComputePlatform.CUDA:
                os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
            else:
                # CPU 模式，不使用 GPU
                if 'CUDA_VISIBLE_DEVICES' in os.environ:
                    del os.environ['CUDA_VISIBLE_DEVICES']

            # 调用 CalvadosWrapper
            wrapper = CalvadosWrapper(temp_config)

            # 写入配置文件到 raw 目录（传入 verbose）
            verbose = self.config.simulation.verbose
            files = wrapper._write_to_dir(raw_dir, gpu_id=gpu_id, verbose=verbose)
            print(f"  📄 配置文件已写入: {files['components']}")
            print(f"  🖥️  平台: {platform.value}")

            # 运行 CALVADOS 模拟
            from multiscale2.extern.ms2_calvados.calvados import sim as calvados_sim
            try:
                calvados_sim.run(
                    path=raw_dir,
                    fconfig='config.yaml',
                    fcomponents='components.yaml'
                )
            except Exception as calvados_error:
                # 如果 CUDA 失败，尝试使用 CPU
                if 'CUDA' in str(calvados_error) or 'Platform' in str(calvados_error) or 'no registered Platform' in str(calvados_error):
                    print(f"  ⚠️  CALVADOS 预平衡 CUDA 失败，尝试使用 CPU...")
                    # 修改配置文件中的 platform 为 CPU
                    import yaml
                    config_file = os.path.join(raw_dir, 'config.yaml')
                    with open(config_file, 'r') as f:
                        config_dict = yaml.safe_load(f)
                    config_dict['platform'] = 'CPU'
                    with open(config_file, 'w') as f:
                        yaml.dump(config_dict, f)
                    
                    # 重新运行
                    calvados_sim.run(
                        path=raw_dir,
                        fconfig='config.yaml',
                        fcomponents='components.yaml'
                    )
                else:
                    raise

            # 查找生成的最终结构
            final_pdb = os.path.join(equilibration_dir, 'final.pdb')
            if os.path.exists(os.path.join(raw_dir, 'checkpoint.pdb')):
                shutil.copy2(
                    os.path.join(raw_dir, 'checkpoint.pdb'),
                    final_pdb
                )
            else:
                # 找带时间戳的 PDB
                for f in os.listdir(raw_dir):
                    if f.endswith('.pdb') and f != 'top.pdb':
                        shutil.copy2(
                            os.path.join(raw_dir, f),
                            final_pdb
                        )
                        break

            # 复制到输出目录根路径（方便后续力场使用）
            output_pdb = os.path.join(self.output_dir, 'preequil_final.pdb')
            if os.path.exists(final_pdb):
                shutil.copy2(final_pdb, output_pdb)
                print(f"  ✓ 预平衡完成: {output_pdb}")
                print(f"  📁 预平衡输出: {equilibration_dir}")
                return output_pdb
            else:
                print(f"  ✗ 未找到预平衡输出文件")
                return None

        except Exception as e:
            print(f"  ✗ 预平衡失败: {e}")
            import traceback
            traceback.print_exc()
            return None

        finally:
            # 恢复原始组件配置
            for i, comp in enumerate(self.config.components):
                orig = original_components[i]
                comp.restraint = orig['restraint']
                comp.restraint_type = orig['restraint_type']
                comp.use_com = orig['use_com']
                comp.k_harmonic = orig['k_harmonic']

    def _create_preequil_config(
        self,
        steps: int = 100000,
        mapping: str = "ca",
        k_restraint: float = 10000.0,
        use_com: bool = True,
        platform: ComputePlatform = ComputePlatform.CUDA,
    ) -> 'CGSimulationConfig':
        """
        创建用于预平衡的临时配置

        对于 MDP 组件：
        - 启用 restraint
        - 设置 restraint_type 为 harmonic
        - 根据 mapping 设置 use_com：
          - CA: use_com=False（约束到每个残基的 CA 原子）
          - COM: use_com=True（约束到每个残基的质心）
        - 设置 k_harmonic 为 k_restraint
        - 步数使用 steps
        - 设置 platform

        Args:
            steps: 预平衡步数
            mapping: 映射方式 ('ca' 或 'com')
            k_restraint: 约束力常数
            use_com: 是否使用 COM 约束
            platform: 计算平台（CUDA 或 CPU）

        Returns:
            临时配置对象
        """
        from copy import deepcopy

        # 深拷贝配置
        temp_config = deepcopy(self.config)

        # 临时修改 MDP 组件的约束设置
        for comp in temp_config.components:
            if comp.type == ComponentType.MDP:
                comp.restraint = True
                comp.restraint_type = 'harmonic'
                comp.use_com = use_com
                comp.k_harmonic = k_restraint

        # 修改模拟参数为预平衡参数
        temp_config.simulation = deepcopy(self.config.simulation)
        temp_config.simulation.steps = steps
        temp_config.simulation.wfreq = min(steps // 10, 1000)
        temp_config.simulation.platform = platform

        return temp_config

    def get_pre_equilibrated_structure(self) -> Optional[str]:
        """
        获取预平衡后的结构文件路径

        Returns:
            预平衡结构文件路径，如果未运行预平衡则返回 None
        """
        if self.output_dir is None:
            return None

        preequil_pdb = os.path.join(self.output_dir, 'preequil_final.pdb')
        if os.path.exists(preequil_pdb):
            return preequil_pdb
        return None

    # -------------------------------------------------------------------------
    # Runner Methods
    # -------------------------------------------------------------------------
    
    def run_calvados(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """
        运行 CALVADOS 模拟

        直接使用 CALVADOS 生成的带时间戳的 PDB 文件作为最终输出。

        Args:
            gpu_id: GPU 设备 ID
            **kwargs: 额外参数

        Returns:
            SimulationResult
        """
        from .calvados_wrapper import CalvadosWrapper
        import time
        from datetime import datetime

        self._ensure_setup()
        self._ensure_not_running()

        # 准备输出目录
        dirs = self.prepare_calvados_output()
        output_dir = dirs['output_dir']
        raw_dir = dirs['raw_dir']
        task_name = dirs['task_name']

        # 复制输入文件（在备份逻辑之后）
        self._copy_input_files(output_dir)

        self.is_running = True
        result = SimulationResult()
        result.output_dir = output_dir

        try:
            print(f"\n[CALVADOS] Running simulation via CGSimulator...")
            print(f"  GPU ID: {gpu_id}")
            print(f"  Task: {task_name}")
            print(f"  Raw output: {raw_dir}")
            print(f"  Topology: {self.config.topol.value if hasattr(self.config.topol, 'value') else self.config.topol}")

            # 设置环境变量指定 GPU
            os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)

            # 调用 CalvadosWrapper 的配置写入和模拟运行（传入 raw_dir、gpu_id 和 verbose）
            verbose = self.config.simulation.verbose
            wrapper = CalvadosWrapper(self.config)
            wrapper._write_to_dir(raw_dir, gpu_id=gpu_id, verbose=verbose)  # 写入配置文件到 raw 目录

            # 运行模拟
            from multiscale2.extern.ms2_calvados.calvados import sim as calvados_sim
            calvados_sim.run(
                path=raw_dir,
                fconfig='config.yaml',
                fcomponents='components.yaml'
            )

            # 复制轨迹文件
            self._copy_trajectory(raw_dir, output_dir)

            # 直接复制带时间戳的 PDB 文件作为 final.pdb
            self._copy_final_pdb(raw_dir, output_dir)

            # 写入日志
            elapsed = 0  # TODO: track actual time
            self._write_simulation_log(output_dir, task_name, elapsed, True)

            result.success = True
            print(f"  ✓ CALVADOS simulation completed")

        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            print(f"  ✗ CALVADOS simulation failed: {e}")

        finally:
            self.is_running = False

        # 设置结果文件路径
        result.trajectory = os.path.join(output_dir, 'trajectory.xtc')
        result.structure = os.path.join(output_dir, 'final.pdb')

        for key in ['trajectory', 'structure']:
            path = getattr(result, key)
            if path and not os.path.exists(path):
                setattr(result, key, None)

        self._result = result
        return result

    def _copy_trajectory(self, raw_dir: str, output_dir: str):
        """从 raw 目录复制轨迹文件"""
        import shutil
        sysname = self.config.system_name
        src_xtc = os.path.join(raw_dir, f'{sysname}.xtc')
        dst_xtc = os.path.join(output_dir, 'trajectory.xtc')
        if os.path.exists(src_xtc):
            shutil.copy2(src_xtc, dst_xtc)
            print(f"  📦 trajectory.xtc")

    def _copy_final_pdb(self, raw_dir: str, output_dir: str):
        """
        从 raw 目录复制带时间戳的 PDB 文件到 output_dir 作为 final.pdb

        CALVADOS 生成的 PDB 文件格式为 {system_name}_{timestamp}.pdb
        """
        import glob

        # 查找带时间戳的 PDB 文件
        pattern = os.path.join(raw_dir, f'{self.config.system_name}_*.pdb')
        pdb_files = glob.glob(pattern)

        if not pdb_files:
            print(f"  ⚠ No timestamped PDB file found in {raw_dir}")
            return

        # 找到最新的文件（按修改时间排序）
        latest_pdb = max(pdb_files, key=os.path.getmtime)

        # 复制到 output_dir/final.pdb
        dst_pdb = os.path.join(output_dir, 'final.pdb')
        shutil.copy2(latest_pdb, dst_pdb)
        print(f"  📦 final.pdb (copied from {os.path.basename(latest_pdb)})")

    def _organize_calvados_output(self, raw_dir: str, output_dir: str, task_name: str):
        """
        整理 CALVADOS 输出文件到统一结构

        统一命名规则：
        - trajectory.xtc  <- {system_name}.xtc
        - final.pdb       <- checkpoint.pdb 或时间戳 PDB
        """
        import shutil

        sysname = self.config.system_name

        # 1. 处理轨迹文件
        src_xtc = os.path.join(raw_dir, f'{sysname}.xtc')
        dst_xtc = os.path.join(output_dir, 'trajectory.xtc')
        if os.path.exists(src_xtc):
            shutil.copy2(src_xtc, dst_xtc)
            print(f"  📦 trajectory.xtc")

        # 2. 查找并复制最终结构
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

        # 3. 复制重要文件
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

    def _write_simulation_log(self, output_dir: str, task_name: str, elapsed: float, success: bool):
        """写入高层级模拟日志"""
        from datetime import datetime

        log_file = os.path.join(output_dir, 'simulation.log')

        status = "SUCCESS" if success else "FAILED"
        components_info = []
        for comp in self.config.components:
            comp_info = f"  - {comp.name}: {comp.type.value if hasattr(comp.type, 'value') else comp.type}, nmol={comp.nmol}"
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
  Topology: {self.config.topol.value if hasattr(self.config.topol, 'value') else self.config.topol}

Components ({len(self.config.components)}):
{chr(10).join(components_info)}

Output Files:
  - final.pdb: Final structure
  - trajectory.xtc: Simulation trajectory
  - raw/: Native simulation output files
"""
        with open(log_file, 'w') as f:
            f.write(log_content)
        print(f"  📝 simulation.log")
    
    def run_hps(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """
        运行 HPS-Urry 模拟

        流程：
        1. 用 CALVADOS 力场构建初始 CG 结构（预平衡）
        2. 使用 HPSParser 解析每个 component（MDP 支持 domain）
        3. 构建 HPSModel 并添加力场
        4. 运行 OpenMM 模拟

        Args:
            gpu_id: GPU 设备 ID
            **kwargs: 额外参数
                - preequil_steps: 预平衡步数
                - preequil_mapping: 预平衡映射方式 ('ca' 或 'com')
                - preequil_k_restraint: 预平衡约束力常数
                - preequil_use_com: 预平衡是否使用 COM 约束
                - preequil_platform: 预平衡平台
                - platform: 模拟平台 (CPU/CUDA)

        Returns:
            SimulationResult
        """
        self._ensure_setup()
        self._ensure_not_running()

        # 预平衡参数
        preequil_steps = kwargs.get('preequil_steps', 100000)
        preequil_mapping = kwargs.get('preequil_mapping', 'ca')
        preequil_k_restraint = kwargs.get('preequil_k_restraint', 10000.0)
        preequil_use_com = kwargs.get('preequil_use_com', True)
        preequil_platform = kwargs.get('preequil_platform', self.config.simulation.platform)
        
        # 检查 CUDA 是否可用，如果不可用则回退到 CPU
        if preequil_platform == ComputePlatform.CUDA:
            try:
                from openmm import Platform
                Platform.getPlatformByName('CUDA')
            except:
                print(f"  ⚠️  CUDA 不可用，预平衡将使用 CPU")
                preequil_platform = ComputePlatform.CPU

        # 预平衡（使用 CALVADOS 构建初始结构）
        preequil_pdb = self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=preequil_steps,
            mapping=preequil_mapping,
            k_restraint=preequil_k_restraint,
            use_com=preequil_use_com,
            platform=preequil_platform,
        )
        if preequil_pdb:
            print(f"  [HPS] 使用预平衡结构: {preequil_pdb}")

        self.is_running = True
        result = SimulationResult()

        try:
            print(f"\n[HPS-Urry] Running simulation...")
            print(f"  GPU ID: {gpu_id}")

            # ===== 1. 导入 HPS 相关模块 =====
            try:
                from multiscale2.extern.ms2_openabc.forcefields.parsers.hps_parser import HPSParser
                from multiscale2.extern.ms2_openabc.forcefields import HPSModel
                from multiscale2.extern.ms2_openabc.lib import _kcal_to_kj
            except ImportError as e:
                raise ImportError(f"ms2_openabc module not available: {e}")

            # ===== 2. 为每个 component 创建 HPSParser =====
            parsers = []
            print("\n  构建 HPSParser...")
            
            for comp in self.config.components:
                if comp.type == ComponentType.MDP:
                    # MDP: 使用提供的 PDB + domain 定义
                    if not comp.fpdb:
                        raise ValueError(f"MDP component '{comp.name}' 需要 fpdb 文件")
                    
                    # 导入原子结构转 CA 的工具函数
                    from multiscale2.extern.ms2_openabc.utils.helper_functions import atomistic_pdb_to_ca_pdb
                    
                    # 创建临时 CA-only PDB 文件（如果原始 PDB 不是 CA-only）
                    ca_pdb_path = comp.fpdb
                    temp_ca_pdb = None
                    
                    # 检查原始 PDB 是否是 CA-only
                    from multiscale2.extern.ms2_openabc.utils import parse_pdb
                    original_atoms = parse_pdb(comp.fpdb)
                    if not (original_atoms['name'].eq('CA').all()):
                        # 原始 PDB 是全原子的，需要转换为 CA-only
                        import os
                        temp_ca_pdb = os.path.join(self.output_dir, f'_temp_{comp.name}_ca.pdb')
                        atomistic_pdb_to_ca_pdb(comp.fpdb, temp_ca_pdb)
                        ca_pdb_path = temp_ca_pdb
                    
                    # 处理 fdomains：可能是文件路径或内联 YAML 内容
                    fdomains_path = comp.fdomains
                    temp_domains_file = None
                    
                    if comp.fdomains:
                        # 检查是否是内联 YAML 内容（包含换行符和 YAML 格式）
                        if '\n' in comp.fdomains and ('[' in comp.fdomains or '-' in comp.fdomains):
                            # 内联 YAML 内容，需要写入临时文件
                            import os
                            import tempfile
                            temp_domains_file = os.path.join(self.output_dir, f'_temp_{comp.name}_domains.yaml')
                            with open(temp_domains_file, 'w') as f:
                                f.write(comp.fdomains)
                            fdomains_path = temp_domains_file
                        elif not os.path.isfile(comp.fdomains):
                            # 既不是文件也不是内联内容，可能是字符串格式的列表
                            try:
                                import ast
                                # 尝试解析为 Python 列表
                                domains_list = ast.literal_eval(comp.fdomains)
                                import os
                                temp_domains_file = os.path.join(self.output_dir, f'_temp_{comp.name}_domains.yaml')
                                # 转换为 YAML 格式
                                import yaml
                                yaml_content = {comp.name: domains_list}
                                with open(temp_domains_file, 'w') as f:
                                    yaml.dump(yaml_content, f)
                                fdomains_path = temp_domains_file
                            except:
                                # 如果解析失败，设为 None 让 HPSParser 处理
                                fdomains_path = None
                    
                    parser = HPSParser(
                        ca_pdb=ca_pdb_path,
                        fdomains=fdomains_path
                    )
                    print(f"    {comp.name}: MDP, {len(parser.atoms)} CA atoms")
                    
                    if parser.enm_pairs:
                        print(f"      → {len(parser.enm_pairs)} ENM pairs")
                    
                    # 清理临时文件
                    if temp_ca_pdb and os.path.exists(temp_ca_pdb):
                        os.remove(temp_ca_pdb)
                    if temp_domains_file and os.path.exists(temp_domains_file):
                        os.remove(temp_domains_file)
                    
                    parsers.append((comp, parser))
                    
                elif comp.type == ComponentType.IDP:
                    # IDP: 从 FASTA 序列构建直的 CA 链
                    if not comp.ffasta:
                        raise ValueError(f"IDP component '{comp.name}' 需要 ffasta 文件")
                    
                    # 读取 FASTA 序列
                    with open(comp.ffasta, 'r') as f:
                        fasta_content = f.read()
                    
                    # 解析 FASTA（提取匹配组件名的序列）
                    lines = fasta_content.strip().split('\n')
                    sequence = None
                    current_seq_lines = []
                    in_target_sequence = False
                    
                    for line in lines:
                        if line.startswith('>'):
                            # 如果上一段是我们要的序列，保存它
                            if in_target_sequence:
                                sequence = ''.join(current_seq_lines)
                                break
                            # 检查这一行是否是我们要找的序列
                            in_target_sequence = (comp.name in line.replace('>', '').strip())
                            current_seq_lines = []
                        elif in_target_sequence:
                            current_seq_lines.append(line.strip())
                    
                    # 处理最后一个序列
                    if sequence is None and in_target_sequence:
                        sequence = ''.join(current_seq_lines)
                    
                    if sequence is None:
                        raise ValueError(f"IDP component '{comp.name}': 未在 FASTA 中找到序列")
                    
                    # 使用 build_straight_CA_chain 从序列构建直的 CA 链
                    from multiscale2.extern.ms2_openabc.utils.helper_functions import build_straight_CA_chain, write_pdb
                    import os
                    
                    # 创建临时 PDB 文件
                    temp_pdb = os.path.join(self.output_dir, f'_temp_idp_{comp.name}_ca.pdb')
                    ca_atoms = build_straight_CA_chain(sequence, r0=0.38)
                    write_pdb(ca_atoms, temp_pdb)
                    
                    # 使用临时 PDB 创建 HPSParser
                    parser = HPSParser(
                        ca_pdb=temp_pdb,
                        fdomains=None
                    )
                    
                    print(f"    {comp.name}: IDP, {len(sequence)} residues (built from FASTA sequence)")
                    parsers.append((comp, parser))
            
            if not parsers:
                raise ValueError("没有有效的 component")

            # ===== 3. 构建 HPSModel =====
            print("\n  构建 HPSModel...")
            model = HPSModel()
            
            # 设置周期性边界
            is_periodic = self.config.topol.value in ['cubic', 'slab']
            model.use_pbc = is_periodic
            
            # 添加所有分子（考虑每个 component 的 nmol）
            for comp, parser in parsers:
                n_before = len(model.atoms) if model.atoms is not None else 0
                # 为每个 component 添加 nmol 个分子
                for _ in range(comp.nmol):
                    model.append_mol(parser)
                n_after = len(model.atoms)
                print(f"    {comp.name}: added {comp.nmol} copies, {n_after - n_before} total atoms")

            # ===== 4. 创建 Topology 和 System =====
            print("\n  创建 Topology 和 System...")
            
            # 创建 HPS 输出子目录（用于存放临时文件和输出文件）
            hps_output_dir = os.path.join(self.output_dir, 'HPS')
            os.makedirs(hps_output_dir, exist_ok=True)
            
            # 从 model.atoms 创建临时 PDB，然后读取创建 topology（放在 HPS 子目录中）
            temp_pdb = os.path.join(hps_output_dir, '_temp_hps_model.pdb')
            model.atoms_to_pdb(temp_pdb, reset_serial=True)
            topology = PDBFile(temp_pdb).topology
            
            # 创建 OpenMM System
            box_size = self.config.box
            is_periodic = self.config.topol.value in ['cubic', 'slab']
            model.create_system(
                top=topology,
                use_pbc=is_periodic,
                box_a=box_size[0],
                box_b=box_size[1],
                box_c=box_size[2]
            )
            print(f"    ✓ System created: {model.system.getNumParticles()} particles")
            
            # ===== 5. 添加力场 =====
            print("\n  添加力场...")
            
            # 蛋白键 (Harmonic)
            print("    - Protein bonds (harmonic)")
            model.add_protein_bonds(force_group=1)
            
            # 非键接触 (Ashbaugh-Hatch with Urry scale)
            print("    - Contacts (Ashbaugh-Hatch, Urry scale)")
            model.add_contacts(
                hydropathy_scale='Urry',
                epsilon=0.2 * _kcal_to_kj,  # Convert kcal to kJ
                mu=1,
                delta=0.08,
                force_group=2
            )
            
            # 静电 (Debye-Hückel)
            print("    - Electrostatics (Debye-Hückel)")
            model.add_dh_elec(
                ldby=1 * unit.nanometer,
                dielectric_water=80.0,
                cutoff=3.5 * unit.nanometer,
                force_group=3
            )
            
            # 弹性网络 (如果有 MDP)
            has_enm = any(p[1].enm_pairs for p in parsers)
            if has_enm:
                print("    - Elastic network (for folded domains)")
                model.add_elastic_network(
                    force_constant=700.0 * unit.kilojoule_per_mole / unit.nanometer ** 2,
                    force_group=4
                )

            # ===== 6. 创建 OpenMM Simulation =====
            print("\n  创建 OpenMM Simulation...")
            
            system = model.system
            
            # 从 CALVADOS 预平衡结构读取初始坐标（与 COCOMO 模式一致）
            if preequil_pdb and os.path.exists(preequil_pdb):
                print(f"  从预平衡结构读取初始坐标: {preequil_pdb}")
                pdb = PDBFile(preequil_pdb)
                positions = pdb.getPositions(asNumpy=True)
                
                # 验证坐标数量是否匹配
                if len(positions) != len(model.atoms):
                    print(f"  ⚠️  警告: 预平衡结构原子数 ({len(positions)}) 与模型原子数 ({len(model.atoms)}) 不匹配")
                    print(f"  从临时 PDB 读取坐标")
                    positions = PDBFile(temp_pdb).getPositions(asNumpy=True)
                else:
                    print(f"  ✓ 坐标数量匹配: {len(positions)} 原子")
            else:
                print(f"  ⚠️  未找到预平衡结构，从临时 PDB 读取坐标")
                positions = PDBFile(temp_pdb).getPositions(asNumpy=True)

            # 选择平台
            config_platform = self.config.simulation.platform
            platform_name = config_platform.value if hasattr(config_platform, 'value') else str(config_platform)
            
            # 尝试 GPU，如果不可用则回退到 CPU
            if gpu_id >= 0:
                for pname in [platform_name, 'CUDA', 'OpenCL']:
                    try:
                        platform = Platform.getPlatformByName(pname)
                        properties = {'DeviceIndex': str(gpu_id)}
                        print(f"    Platform: {pname}")
                        break
                    except:
                        continue
                else:
                    platform = Platform.getPlatformByName('CPU')
                    properties = {}
                    print(f"    Platform: CPU (GPU unavailable)")
            else:
                platform = Platform.getPlatformByName('CPU')
                properties = {}
                print(f"    Platform: CPU")

            # 创建积分器
            temperature = self.config.temperature
            integrator = LangevinMiddleIntegrator(
                temperature * kelvin,
                0.01 / picoseconds,
                0.01 * picoseconds
            )

            simulation = Simulation(
                topology,
                system,
                integrator=integrator,
                platform=platform,
                platformProperties=properties
            )

            # 设置初始坐标
            simulation.context.setPositions(positions)

            # 设置盒子向量
            if is_periodic:
                box_size = self.config.box
                box_vecs = [
                    mm.Vec3(x=box_size[0], y=0.0, z=0.0),
                    mm.Vec3(x=0.0, y=box_size[1], z=0.0),
                    mm.Vec3(x=0.0, y=0.0, z=box_size[2])
                ] * unit.nanometer
                simulation.context.setPeriodicBoxVectors(*box_vecs)

            # 能量最小化
            print("  Energy minimization...")
            simulation.minimizeEnergy()
            print("  Energy minimization completed")
            simulation.context.setVelocitiesToTemperature(temperature * kelvin)

            # ===== 6. 运行模拟 =====
            print(f"\n  Running HPS simulation: {self.config.simulation.steps} steps...")
            
            # 添加报告器（使用之前创建的 HPS 输出目录）
            output_dir = hps_output_dir
            
            wfreq = self.config.simulation.wfreq
            traj_file = os.path.join(output_dir, 'trajectory.xtc')
            log_file = os.path.join(output_dir, 'simulation.log')
            
            simulation.reporters.append(XTCReporter(traj_file, wfreq))
            simulation.reporters.append(StateDataReporter(
                log_file,
                wfreq,
                step=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                volume=True,
            ))

            # 进度条
            from tqdm import tqdm
            total_steps = self.config.simulation.steps
            batch_size = 1000
            
            for _ in tqdm(range(0, total_steps, batch_size), desc="HPS"):
                simulation.step(min(batch_size, total_steps - simulation.currentStep))

            print("  Simulation completed!")

            # ===== 7. 保存结果 =====
            print("\n  保存结果...")
            
            # 最终结构
            state_final = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=True,
                getEnergy=True,
                enforcePeriodicBox=True
            )
            positions_final = state_final.getPositions()
            
            final_pdb = os.path.join(output_dir, 'final.pdb')
            with open(final_pdb, 'w') as f:
                PDBFile.writeFile(topology, positions_final, f)
            print(f"    - final.pdb")
            
            # 复制到根目录（output_dir 是 HPS/，需要复制到 TDP43_CTD_CG/）
            system_name = self.config.system_name
            cg_root_dir = os.path.dirname(output_dir)  # output_dir 是 .../TDP43_CTD_CG/HPS，所以 dirname 是 TDP43_CTD_CG
            final_pdb_root = os.path.join(cg_root_dir, 'final.pdb')
            os.makedirs(cg_root_dir, exist_ok=True)
            shutil.copy2(final_pdb, final_pdb_root)
            print(f"  复制最终结构到: {final_pdb_root}")

            # 保存系统 XML
            system_xml = os.path.join(output_dir, 'system.xml')
            with open(system_xml, 'w') as f:
                f.write(XmlSerializer.serialize(system))
            print(f"    - system.xml")

            # 保存 checkpoint
            checkpoint_file = os.path.join(output_dir, 'restart.chk')
            simulation.saveCheckpoint(checkpoint_file)
            print(f"    - restart.chk")

            result.success = True
            result.trajectory = traj_file
            result.structure = final_pdb
            result.output_dir = output_dir
            
            print(f"\n  HPS-Urry 输出目录: {output_dir}")

        except ImportError as e:
            result.success = False
            result.errors.append(f"ms2_openabc module not available: {e}")
            print(f"  ✗ HPS-Urry simulation failed: ms2_openabc module not available")
            import traceback
            traceback.print_exc()
        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            print(f"  ✗ HPS-Urry simulation failed: {e}")
            import traceback
            traceback.print_exc()

        finally:
            self.is_running = False

        self._result = result
        return result

    def run_moff(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """
        运行 MOFF 模拟

        自动进行预平衡（使用 CALVADOS 构建初始结构）：
        1. 用 CALVADOS 力场构建初始 CG 结构
        2. 对 MDP 蛋白进行短时间模拟
        3. 然后切换到 MOFF 力场进行正式模拟

        使用 OpenABC 包的 MOFF 力场。

        预平衡参数（硬编码，可通过 kwargs 覆盖）：
        - steps: 预平衡步数（默认 100000）
        - mapping: 映射方式 'ca' 或 'com'（默认 'ca'）
        - k_restraint: 约束力常数（默认 10000.0）
        - use_com: 是否使用 COM 约束（默认 True）
        - platform: 计算平台（默认从 config 读取，CUDA）

        Args:
            gpu_id: GPU 设备 ID
            **kwargs: 额外参数
                - salt_conc: 盐浓度（默认 150 mM）
                - preequil_steps: 预平衡步数
                - preequil_mapping: 预平衡映射方式 ('ca' 或 'com')
                - preequil_k_restraint: 预平衡约束力常数
                - preequil_use_com: 预平衡是否使用 COM 约束
                - preequil_platform: 预平衡平台（CUDA 或 CPU）

        Returns:
            SimulationResult
        """
        self._ensure_setup()
        self._ensure_not_running()

        # 预平衡参数（硬编码，可覆盖）
        preequil_steps = kwargs.get('preequil_steps', 100000)
        preequil_mapping = kwargs.get('preequil_mapping', 'ca')
        preequil_k_restraint = kwargs.get('preequil_k_restraint', 10000.0)
        preequil_use_com = kwargs.get('preequil_use_com', True)
        # 预平衡平台：优先使用 kwargs 中指定的值，否则使用 config 中的设置（默认为 CUDA）
        preequil_platform = kwargs.get('preequil_platform', self.config.simulation.platform)

        # 预平衡（使用 CALVADOS 构建初始结构）
        preequil_pdb = self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=preequil_steps,
            mapping=preequil_mapping,
            k_restraint=preequil_k_restraint,
            use_com=preequil_use_com,
            platform=preequil_platform,
        )
        if preequil_pdb:
            print(f"  [MOFF] 使用预平衡结构: {preequil_pdb}")

        self.is_running = True
        result = SimulationResult()
        result.output_dir = self.output_dir

        try:
            print(f"\n[MOFF] Running simulation...")
            print(f"  GPU ID: {gpu_id}")

            # TODO: 实现 MOFF runner
            # 需要使用 openabc.forcefields.MOFFModel

            result.success = True
            print(f"  ✓ MOFF simulation completed (placeholder)")

        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            print(f"  ✗ MOFF simulation failed: {e}")

        finally:
            self.is_running = False

        self._result = result
        return result
    
    def run_cocomo(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """
        运行 COCOMO 模拟

        自动进行预平衡（使用 CALVADOS 构建初始结构）：
        1. 用 CALVADOS 力场构建初始 CG 结构
        2. 对 MDP 蛋白进行短时间模拟
        3. 然后切换到 COCOMO 力场进行正式模拟

        COCOMO 使用 COCOMO2 力场，默认启用 SASA 修正。
        离子强度固定为 0.1 M（力场默认值）。

        预平衡参数（硬编码，可通过 kwargs 覆盖）：
        - steps: 预平衡步数（默认 100000）
        - mapping: 映射方式 'ca' 或 'com'（默认 'ca'）
        - k_restraint: 约束力常数（默认 10000.0）
        - use_com: 是否使用 COM 约束（默认 True）
        - platform: 计算平台（默认从 config 读取，CUDA）

        Args:
            gpu_id: GPU 设备 ID
            **kwargs: 额外参数
                - nstep: 模拟步数（覆盖配置中的值）
                - wfreq: 写入频率
                - tstep: 时间步长（ps，默认 0.01）
                - gamma: 摩擦系数（1/ps，默认 0.01）
                - preequil_steps: 预平衡步数
                - preequil_mapping: 预平衡映射方式 ('ca' 或 'com')
                - preequil_k_restraint: 预平衡约束力常数
                - preequil_use_com: 预平衡是否使用 COM 约束
                - preequil_platform: 预平衡平台（CUDA 或 CPU）

        Returns:
            SimulationResult
        """
        import warnings as py_warnings

        self._ensure_setup()
        self._ensure_not_running()

        # Warning for ionic strength
        if abs(self.config.ionic - 0.1) > 0.001:
            py_warnings.warn(
                "COCOMO uses default ionic strength of 0.1 M. "
                f"Config value {self.config.ionic} M will be ignored.",
                UserWarning
            )

        # 预平衡参数（硬编码，可覆盖）
        preequil_steps = kwargs.get('preequil_steps', 100000)
        preequil_mapping = kwargs.get('preequil_mapping', 'ca')
        preequil_k_restraint = kwargs.get('preequil_k_restraint', 10000.0)
        preequil_use_com = kwargs.get('preequil_use_com', True)
        # 预平衡平台：优先使用 kwargs 中指定的值，否则使用 config 中的设置（默认为 CUDA）
        preequil_platform = kwargs.get('preequil_platform', self.config.simulation.platform)

        # 预平衡（使用 CALVADOS 构建初始结构）
        preequil_pdb = self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=preequil_steps,
            mapping=preequil_mapping,
            k_restraint=preequil_k_restraint,
            use_com=preequil_use_com,
            platform=preequil_platform,
        )

        # 准备输出目录
        dirs = self._prepare_cocomo_output()
        output_dir = dirs['output_dir']

        self.is_running = True
        result = SimulationResult()

        try:
            print(f"\n[COCOMO] Running COCOMO2 simulation...")
            print(f"  预平衡结构: {preequil_pdb}")

            # 使用新的 COCOMO 类运行模拟
            from .cocomo2_creator import COCOMO

            # 构建 topology_info 字典
            topology_info = {
                'global_sequence': self.get_global_sequence(),
                'chain_ids': self.get_chain_ids(),
                'folded_domains': self.get_folded_domains(),
                'component_names': self._build_component_names(),
                'local_residue_indices': list(range(1, len(self.get_global_sequence()) + 1)),
            }

            # 读取预平衡结构的坐标
            from openmm.app import PDBFile
            pdb = PDBFile(preequil_pdb)
            positions = pdb.getPositions(asNumpy=True)

            # 获取 box size（使用 config 中的 box 参数，单位为 nm）
            box_size = self.config.box

            # 获取 SASA 值
            sasa_values = self._get_sasa_values()

            if sasa_values is not None:
                topology_info['sasa_values'] = sasa_values

            # 创建 COCOMO 系统
            cocomo = COCOMO(
                box_size=box_size,
                topology_info=topology_info,
                positions=positions,
                surf=0.7,
                resources='CUDA' if gpu_id >= 0 else 'CPU'
            )

            # 设置 _topology_info 以便 ENM 计算使用
            # 使用 SimpleNamespace 创建带有属性的对象
            from types import SimpleNamespace
            cocomo._topology_info = SimpleNamespace(
                global_sequence=topology_info['global_sequence'],
                chain_ids=topology_info['chain_ids'],
                is_folded=topology_info['folded_domains'],
                sasa_values=topology_info.get('sasa_values', [])
            )

            # 创建 OpenMM 系统
            system, topology = cocomo.create_system()

            # 创建 Simulation 对象 (优先使用config指定的platform)
            config_platform = self.config.simulation.platform
            platform_name = config_platform.value if hasattr(config_platform, 'value') else str(config_platform)
            
            # 平台选择：优先使用config指定的平台，支持自动回退
            if gpu_id >= 0:
                # 用户想要使用GPU
                try:
                    platform = Platform.getPlatformByName(platform_name)
                    properties = {'DeviceIndex': str(gpu_id)}
                    print(f"  使用 {platform_name} 平台")
                except Exception:
                    # 尝试其他GPU平台
                    for fallback in ['CUDA', 'OpenCL']:
                        if fallback != platform_name:
                            try:
                                platform = Platform.getPlatformByName(fallback)
                                properties = {'DeviceIndex': str(gpu_id)}
                                print(f"  {platform_name} 不可用，回退到 {fallback}")
                                break
                            except:
                                continue
                    else:
                        print(f"  GPU 不可用，回退到 CPU")
                        platform = Platform.getPlatformByName('CPU')
                        properties = {}
            else:
                # 用户想要使用CPU
                platform = Platform.getPlatformByName('CPU')
                properties = {}
                print("  使用 CPU 平台")

            # 使用 LangevinIntegrator（与 Legacy 版本一致）
            # 温度从 config 读取，摩擦系数 0.01/ps，时间步长 0.01 ps
            temperature = self.config.temperature
            integrator = LangevinIntegrator(
                temperature * kelvin,
                0.01 / picoseconds,
                0.01 * picoseconds
            )

            simulation = Simulation(
                topology,
                system,
                integrator=integrator,
                platform=platform,
                platformProperties=properties
            )

            # 设置初始坐标
            simulation.context.setPositions(positions)

            # 设置初始速度（与 Legacy 版本一致）
            simulation.context.setVelocitiesToTemperature(temperature * kelvin)

            # 能量最小化
            print("  Running energy minimization...")
            simulation.minimizeEnergy()
            print("  Energy minimization completed")

            # 运行模拟
            wfreq = self.config.simulation.wfreq
            xtc_file = os.path.join(output_dir, 'trajectory.xtc')
            log_file = os.path.join(output_dir, 'simulation.log')

            # 添加报告器（使用 mdtraj 的 XTCReporter 保留 PBC 信息）
            simulation.reporters.append(
                XTCReporter(xtc_file, wfreq)
            )
            simulation.reporters.append(
                StateDataReporter(
                    log_file,
                    wfreq,
                    step=True,
                    potentialEnergy=True,
                    kineticEnergy=True,
                    totalEnergy=True,
                    temperature=True,
                    volume=True,
                )
            )

            # 运行模拟（使用 tqdm 显示进度）
            total_steps = self.config.simulation.steps
            print(f"  开始模拟: {total_steps} 步")

            from tqdm import tqdm
            n_batches = 10
            batch_size = total_steps // n_batches

            for _ in tqdm(range(n_batches), desc="COCOMO"):
                simulation.step(batch_size)
                simulation.saveCheckpoint(os.path.join(output_dir, 'restart.chk'))

            # 处理剩余步数
            remaining = total_steps % n_batches
            if remaining > 0:
                simulation.step(remaining)

            print(f"  模拟完成!")

            # 获取最终状态（包含 PBC 信息）
            state_final = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=True,
                getEnergy=True,
                enforcePeriodicBox=True  # 强制周期边界条件
            )

            # 获取最终的位置和盒子向量
            positions_final = state_final.getPositions()
            box_vectors = state_final.getPeriodicBoxVectors()

            # 在拓扑上设置盒子向量（这样 PDB 会包含 PBC 信息）
            simulation.topology.setPeriodicBoxVectors(box_vectors)

            # 保存最终结构为 PDB 格式（包含 PBC）
            from openmm.app import PDBFile
            final_pdb = os.path.join(output_dir, 'final.pdb')
            with open(final_pdb, 'w') as f:
                PDBFile.writeFile(
                    simulation.topology,
                    positions_final,
                    f
                )
            print(f"  保存最终结构: {final_pdb}")

            # final.pdb 已经在根目录（output_dir = self.output_dir = {system_name}_CG），无需复制

            # 保存系统 XML
            system_xml = os.path.join(output_dir, 'system.xml')
            with open(system_xml, 'w') as f:
                f.write(XmlSerializer.serialize(system))
            print(f"  保存系统 XML: {system_xml}")

            result.success = True
            result.trajectory = xtc_file
            result.structure = final_pdb
            result.output_dir = output_dir

            print(f"  COCOMO 输出目录: {output_dir}")

        except Exception as e:
            result = SimulationResult()
            result.success = False
            result.errors.append(str(e))
            result.output_dir = output_dir
            print(f"  ✗ COCOMO simulation failed: {e}")
            import traceback
            traceback.print_exc()

        finally:
            self.is_running = False

        self._result = result
        return result
    
    def _prepare_cocomo_output(self) -> Dict[str, str]:
        """
        准备 COCOMO 输出目录

        输出目录结构：{output_dir}/ (CLI already sets self.output_dir = {system_name}_CG)
        直接使用 self.output_dir，避免嵌套目录
        """
        # CLI 已经设置了 self.output_dir = {system_name}_CG，直接使用
        output_dir = self.output_dir

        # 备份旧结果（仅在目录包含之前的模拟结果时）
        import shutil
        from datetime import datetime

        # 检查目录是否存在，以及是否包含之前的模拟结果
        should_backup = False
        if os.path.exists(output_dir):
            # 检查是否有模拟结果文件（不是预平衡文件）
            result_files = ['final.pdb', 'trajectory.xtc', 'simulation.log', 'system.xml']
            has_results = any(os.path.exists(os.path.join(output_dir, f)) for f in result_files)
            
            # 只有当存在模拟结果时才备份（预平衡文件不算）
            if has_results:
                should_backup = True
        
        if should_backup:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{output_dir}_backup_{timestamp}"
            shutil.move(output_dir, backup_dir)
            print(f"  📁 备份旧结果到: {backup_dir}")
            # 重新创建输出目录
            os.makedirs(output_dir, exist_ok=True)
        else:
            # 目录不存在或只包含预平衡文件，直接创建（如果不存在）
            os.makedirs(output_dir, exist_ok=True)

        return {
            'output_dir': output_dir,
            'task_name': 'COCOMO',
        }

    def run_mpipi_recharged(self, gpu_id: int = 0, **kwargs) -> SimulationResult:
        """
        运行 Mpipi-Recharged 模拟

        使用与 test_100_molecules.py 相同的逻辑：
        1. 创建 biomolecule 对象（MDP/IDP）
        2. 使用 gmx insert-molecules 放置分子
        3. 构建 Mpipi-Recharged 系统
        4. 能量最小化
        5. 运行 MD 模拟

        注意：此方法不使用 CALVADOS 预平衡，Mpipi-Recharged 自己会处理所有初始化。

        Args:
            gpu_id: GPU 设备 ID
            **kwargs: 额外参数
                - use_gmx_insert: 是否使用 gmx insert-molecules 放置分子（默认 True）
                - gmx_radius: gmx insert-molecules 的最小距离（nm，默认 0.35）
                - verbose: 是否输出详细信息（默认 True）

        Returns:
            SimulationResult
        """
        from openmm.app import PDBFile
        from openmm import Platform, LangevinMiddleIntegrator, Vec3, LocalEnergyMinimizer
        import openmm as mm
        from openmm import unit
        from openmm.unit import kelvin, picoseconds, nanometer

        self._ensure_setup()
        self._ensure_not_running()

        # 输出目录
        dirs = self._prepare_mpipi_output()
        output_dir = dirs['output_dir']

        self.is_running = True
        result = SimulationResult()

        try:
            print(f"\n[Mpipi-Recharged] Running simulation...")

            # 从 ms2_OpenMpipi 导入 biomolecule 类和新函数
            from multiscale2.extern.ms2_OpenMpipi import MDP, IDP, build_mpipi_recharged_system_from_chains

            # 从 config.components 构建 biomolecule 对象
            print("\n  构建 biomolecule 对象...")
            chain_objects = []
            for comp in self.config.components:
                # 获取序列
                if comp.type == ComponentType.IDP:
                    # IDP 从 fasta 读取序列
                    if comp.ffasta:
                        with open(comp.ffasta, 'r') as f:
                            fasta_content = f.read()
                        # 解析 FASTA（只提取匹配组件名的序列）
                        lines = fasta_content.strip().split('\n')
                        sequence = None
                        current_seq_lines = []
                        in_target_sequence = False
                        
                        for line in lines:
                            if line.startswith('>'):
                                # 如果上一段是我们要的序列，保存它
                                if in_target_sequence:
                                    sequence = ''.join(current_seq_lines)
                                    break
                                # 检查这一行是否是我们要找的序列
                                in_target_sequence = (comp.name in line.replace('>', '').strip())
                                current_seq_lines = []
                            elif in_target_sequence:
                                current_seq_lines.append(line.strip())
                        
                        # 处理最后一个序列
                        if sequence is None and in_target_sequence:
                            sequence = ''.join(current_seq_lines)
                        
                        if sequence is None:
                            print(f"    ⚠️  {comp.name}: 未在 FASTA 中找到序列，跳过")
                            continue
                    else:
                        continue
                    
                    # 创建 IDP 对象
                    idp = IDP(comp.name, sequence)
                    chain_objects.append((comp, idp))
                    print(f"    {comp.name}: IDP, {len(sequence)} residues")

                elif comp.type == ComponentType.MDP:
                    # MDP 从 fpdb 读取结构和序列
                    if not comp.fpdb:
                        continue
                    
                    # 解析 fdomains 获取折叠域索引
                    domains = self._parse_fdomains(comp.fdomains)
                    
                    # 从 PDB 读取序列（MDP 永远使用 PDB 中的序列）
                    # 注意：PDB 残基名是三字母码，需要转换为单字母码
                    pdb_temp = PDBFile(comp.fpdb)
                    three_to_one = {
                        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                    }
                    sequence = []
                    for res in pdb_temp.topology.residues():
                        if res.name in three_to_one:
                            sequence.append(three_to_one[res.name])
                    sequence = ''.join(sequence)
                    
                    # 创建 MDP 对象
                    mdp = MDP(comp.name, sequence, domains, comp.fpdb)
                    chain_objects.append((comp, mdp))
                    print(f"    {comp.name}: MDP, {len(sequence)} residues, {len(domains)} domains")

            if not chain_objects:
                raise ValueError("No valid biomolecules could be constructed for Mpipi simulation.")

            # 构建 chain_info 字典
            print("\n  构建 chain_info 字典...")
            chain_info = {}
            for comp, biomol in chain_objects:
                chain_info[biomol] = comp.nmol
                print(f"    {comp.name}: {comp.nmol} copies")

            # 调用 build_mpipi_recharged_system_from_chains
            # 这个函数会：1) relax 每个单体 2) build model 3) 添加力场
            csx = self.config.ionic * 1000  # M -> mM
            is_periodic = self.config.topol.value in ['cubic', 'slab']
            box_size = self.config.box  # 使用 config 指定的盒子大小
            
            print("\n  构建 Mpipi-Recharged 系统（包含 relaxation + model building + 力场）...")
            # 使用 gmx insert-molecules 作为默认的分子放置方法（更稳健，避免重叠）
            use_gmx_insert = kwargs.get('use_gmx_insert', True)  # 默认使用 gmx insert-molecules
            gmx_radius = kwargs.get('gmx_radius', 0.35)  # 默认最小距离 0.35 nm
            verbose = kwargs.get('verbose', True)
            
            if use_gmx_insert:
                print(f"  使用 gmx insert-molecules 放置分子 (radius={gmx_radius} nm)")
            else:
                print(f"  使用 CALVADOS 网格放置方法")
            
            system, model = build_mpipi_recharged_system_from_chains(
                chain_info=chain_info,
                box_size=box_size,
                topol=self.config.topol.value,
                T=self.config.temperature * unit.kelvin,
                csx=csx,
                CM_remover=True,
                periodic=is_periodic,
                use_gmx_insert=use_gmx_insert,  # 使用 gmx insert-molecules
                gmx_radius=gmx_radius,  # 最小距离
                verbose=verbose
            )
            print(f"  系统构建完成: {system.getNumParticles()} 粒子, {system.getNumForces()} 力")

            # 从 model 获取 topology 和 positions
            topology = model.topology
            positions = model.positions

            # 创建 Simulation 对象 (优先使用config指定的platform)
            config_platform = self.config.simulation.platform
            platform_name = config_platform.value if hasattr(config_platform, 'value') else str(config_platform)
            
            # 平台选择：优先使用config指定的平台，支持自动回退
            if gpu_id >= 0:
                # 用户想要使用GPU
                try:
                    platform = Platform.getPlatformByName(platform_name)
                    properties = {'DeviceIndex': str(gpu_id)}
                    print(f"  使用 {platform_name} 平台 (GPU {gpu_id})")
                except Exception:
                    # 尝试其他GPU平台
                    for fallback in ['CUDA', 'OpenCL']:
                        if fallback != platform_name:
                            try:
                                platform = Platform.getPlatformByName(fallback)
                                properties = {'DeviceIndex': str(gpu_id)}
                                print(f"  {platform_name} 不可用，回退到 {fallback}")
                                break
                            except:
                                continue
                    else:
                        print(f"  GPU 不可用，回退到 CPU")
                        platform = Platform.getPlatformByName('CPU')
                        properties = {}
            else:
                # 用户想要使用CPU
                platform = Platform.getPlatformByName('CPU')
                properties = {}
                print("  使用 CPU 平台")

            # 使用 LangevinMiddleIntegrator（与 test_100_molecules.py 一致）
            temperature = self.config.temperature
            integrator = LangevinMiddleIntegrator(
                temperature * kelvin,
                0.01 / picoseconds,
                0.01 * picoseconds
            )

            simulation = mm.app.Simulation(
                topology,
                system,
                integrator=integrator,
                platform=platform,
                platformProperties=properties
            )

            # 设置初始坐标：使用 gmx insert-molecules 生成的坐标（不使用 CALVADOS 预平衡）
            print(f"  使用 gmx insert-molecules 生成的坐标")
            simulation.context.setPositions(positions)

            # 如果系统有周期性边界，设置盒子向量
            if is_periodic:
                box_size = self.config.box
                box_vecs = [
                    mm.Vec3(x=box_size[0], y=0.0, z=0.0),
                    mm.Vec3(x=0.0, y=box_size[1], z=0.0),
                    mm.Vec3(x=0.0, y=0.0, z=box_size[2])
                ] * unit.nanometer
                simulation.context.setPeriodicBoxVectors(*box_vecs)

            # 保存 minimization 前的结构（用于调试）
            debug_pdb_before_min = os.path.join(output_dir, 'before_minimization.pdb')
            with open(debug_pdb_before_min, 'w') as f:
                PDBFile.writeFile(
                    simulation.topology,
                    positions,
                    f
                )
            print(f"  保存 minimization 前结构: {debug_pdb_before_min}")

            # 能量最小化（与 test_100_molecules.py 一致：minimizeEnergy 不带参数）
            print("  Running energy minimization...")
            simulation.minimizeEnergy()
            print("  Energy minimization completed")

            # 获取最小化后的状态和能量
            state_after_min = simulation.context.getState(getPositions=True, getEnergy=True)
            energy_after_min = state_after_min.getPotentialEnergy()
            print(f"  Minimization energy: {energy_after_min}")

            # 检查能量是否为 NaN
            if np.isnan(energy_after_min.value_in_unit(unit.kilojoule_per_mole)):
                raise RuntimeError("Energy minimization failed: NaN energy")

            # 保存最小化后的结构
            positions_after_min = state_after_min.getPositions()
            debug_pdb_after_min = os.path.join(output_dir, 'after_minimization.pdb')
            with open(debug_pdb_after_min, 'w') as f:
                PDBFile.writeFile(
                    simulation.topology,
                    positions_after_min,
                    f
                )
            print(f"  保存 minimization 后结构: {debug_pdb_after_min}")

            # 运行模拟
            wfreq = self.config.simulation.wfreq
            xtc_file = os.path.join(output_dir, 'trajectory.xtc')
            log_file = os.path.join(output_dir, 'simulation.log')

            # 添加报告器
            simulation.reporters.append(
                XTCReporter(xtc_file, wfreq)
            )
            simulation.reporters.append(
                StateDataReporter(
                    log_file,
                    wfreq,
                    step=True,
                    potentialEnergy=True,
                    kineticEnergy=True,
                    totalEnergy=True,
                    temperature=True,
                    volume=True,
                )
            )

            # 运行模拟（使用 tqdm 显示进度）
            total_steps = self.config.simulation.steps
            print(f"  开始模拟: {total_steps} 步")

            from tqdm import tqdm
            n_batches = 10
            batch_size = total_steps // n_batches

            for _ in tqdm(range(n_batches), desc="Mpipi-Recharged"):
                simulation.step(batch_size)
                simulation.saveCheckpoint(os.path.join(output_dir, 'restart.chk'))

            # 处理剩余步数
            remaining = total_steps % n_batches
            if remaining > 0:
                simulation.step(remaining)

            print(f"  模拟完成!")

            # 获取最终状态（包含 PBC 信息）
            state_final = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=True,
                getEnergy=True,
                enforcePeriodicBox=True
            )

            # 获取最终的位置和盒子向量
            positions_final = state_final.getPositions()
            box_vectors = state_final.getPeriodicBoxVectors()

            # 在拓扑上设置盒子向量
            simulation.topology.setPeriodicBoxVectors(box_vectors)

            # 保存最终结构为 PDB 格式（包含成键信息 - CONECT 记录）
            final_pdb = os.path.join(output_dir, 'final.pdb')
            with open(final_pdb, 'w') as f:
                PDBFile.writeFile(
                    simulation.topology,
                    positions_final,
                    f,
                    keepIds=True
                )

            # 复制 final.pdb 到根目录（self.output_dir = {system_name}_CG）
            import shutil
            final_pdb_root = os.path.join(self.output_dir, 'final.pdb')
            shutil.copy2(final_pdb, final_pdb_root)
            print(f"  复制最终结构到根目录: {final_pdb_root}")

            # 复制文件到标准位置
            shutil.copy(final_pdb, os.path.join(output_dir, 'preequil_final.pdb'))

            # 保存系统 XML
            system_xml = os.path.join(output_dir, 'system.xml')
            with open(system_xml, 'w') as f:
                f.write(mm.XmlSerializer.serialize(system))

            # 保存检查点
            shutil.copy(os.path.join(output_dir, 'restart.chk'), os.path.join(output_dir, 'final.chk'))

            # 打印结果摘要
            print(f"\n  ✓ Simulation completed!")
            print(f"  Final structure: {final_pdb}")
            print(f"  Trajectory: {xtc_file}")
            
            # 打印最终能量
            final_energy = state_final.getPotentialEnergy()
            print(f"  Final potential energy: {final_energy}")
            print(f"  Final energy per particle: {final_energy.value_in_unit(unit.kilojoule_per_mole) / system.getNumParticles():.3f} kJ/mol")

            result.success = True
            result.output_dir = output_dir
            result.trajectory = xtc_file
            result.structure = final_pdb
            result.errors = []

        except Exception as e:
            import traceback
            print(f"\n  ✗ Mpipi-Recharged simulation failed: {e}")
            if kwargs.get('verbose', False):
                traceback.print_exc()
            
            result.success = False
            result.output_dir = output_dir
            result.errors = [str(e)]

        finally:
            self.is_running = False

        return result
    
    def _build_globular_indices_dict(self) -> Dict[str, list]:
        """
        从 config.components 构建 globular_indices_dict

        用于 OpenMpipi 的 get_mpipi_system 函数。
        OpenMpipi 期望格式：{chain_id: [[start1, end1], [start2, end2], ...]}
        其中每个 [start, end] 是一个折叠域的索引范围（inclusive）。
        
        注意：OpenMpipi 期望的是**局部索引**（相对于每个chain的0-based索引），
        而不是全局系统索引。

        Returns:
            Dictionary mapping chain_id to list of domain ranges [start, end] (local indices)
        """
        globular_indices_dict = {}

        # 获取链ID列表
        chain_ids = self.get_chain_identifiers()

        # 获取 folded domain 信息
        folded_domains = self.get_folded_domains()

        # 构建字典：{chain_id: [[start1, end1], [start2, end2], ...]}
        # 首先收集所有链的起始位置（局部索引的基准）
        chain_local_start = {}  # {chain_id: local_index_offset}
        
        for res_idx, chain_id in enumerate(chain_ids):
            if chain_id not in chain_local_start:
                chain_local_start[chain_id] = res_idx

        # 遍历每个残基，检测连续的 folded 区域（域范围）
        current_chain = None
        domain_start = None  # 全局索引
        domain_start_local = None  # 局部索引

        for res_idx, (chain_id, is_folded) in enumerate(zip(chain_ids, folded_domains)):
            if chain_id not in globular_indices_dict:
                globular_indices_dict[chain_id] = []

            # 新链开始，重置状态
            if current_chain != chain_id:
                if current_chain is not None and domain_start is not None:
                    # 保存上一个域（使用局部索引）
                    local_start = domain_start - chain_local_start[current_chain]
                    local_end = (res_idx - 1) - chain_local_start[current_chain]
                    globular_indices_dict[current_chain].append([local_start, local_end])
                current_chain = chain_id
                domain_start = None

            # 如果是 folded domain，记录起始位置
            if is_folded:
                if domain_start is None:
                    domain_start = res_idx
            else:
                # 如果之前在域中，现在结束了，保存域范围
                if domain_start is not None:
                    local_start = domain_start - chain_local_start[chain_id]
                    local_end = (res_idx - 1) - chain_local_start[chain_id]
                    globular_indices_dict[chain_id].append([local_start, local_end])
                    domain_start = None

        # 处理最后一个域（如果链末尾是 folded）
        if current_chain is not None and domain_start is not None:
            local_start = domain_start - chain_local_start[current_chain]
            # 局部索引的结束位置是该链的最后一个残基
            chain_length = sum(1 for cid in chain_ids if cid == current_chain)
            local_end = chain_length - 1
            globular_indices_dict[current_chain].append([local_start, local_end])

        return globular_indices_dict

    def _prepare_mpipi_output(self) -> Dict[str, str]:
        """
        准备 Mpipi-Recharged 输出目录

        输出目录结构：{output_dir}/Mpipi-Recharged/
        CLI 已经设置了 self.output_dir = {system_name}_CG，直接使用

        期望结构：
        {system_name}_CG/
        ├── Mpipi-Recharged/     # 主模拟输出
        │   ├── trajectory.xtc
        │   ├── final.pdb
        │   └── ...
        ├── final.pdb            # 复制到根目录
        └── equilibration/       # 预平衡输出（由 _run_pre_equilibration 创建）
            └── raw/
                └── ...
        """
        # CLI 已经设置了 self.output_dir = {system_name}_CG，直接使用
        # Mpipi-Recharged 主模拟输出目录
        mpipi_dir = os.path.join(self.output_dir, 'Mpipi-Recharged')

        # 备份旧结果
        import shutil
        from datetime import datetime

        if os.path.exists(mpipi_dir):
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{mpipi_dir}_backup_{timestamp}"
            shutil.move(mpipi_dir, backup_dir)
            print(f"  📁 备份旧结果到: {backup_dir}")

        os.makedirs(mpipi_dir, exist_ok=True)

        return {
            'output_dir': mpipi_dir,
            'task_name': 'Mpipi-Recharged',
        }
    
    # -------------------------------------------------------------------------
    # Utility Methods
    # -------------------------------------------------------------------------

    def get_result(self) -> Optional[SimulationResult]:
        """获取最近的模拟结果"""
        return self._result

    def cleanup(self):
        """清理临时文件"""
        self.is_setup = False
        self._result = None

    # -------------------------------------------------------------------------
    # Pre-equilibration Utility Methods
    # -------------------------------------------------------------------------

    def run_pre_equilibration(
        self,
        gpu_id: int = 0,
        steps: int = 100000,
        mapping: str = "ca",
        k_restraint: float = 10000.0,
        use_com: bool = True,
        platform: Optional[ComputePlatform] = None,
    ) -> Optional[str]:
        """
        单独运行预平衡（使用 CALVADOS 构建初始结构）

        此方法允许用户在不运行完整模拟的情况下，预先生成 CG 结构。
        生成的 `preequil_final.pdb` 文件可用于后续的力场模拟。

        Args:
            gpu_id: GPU 设备 ID
            steps: 预平衡步数（默认 100000）
            mapping: 映射方式 ('ca' 或 'com')（默认 'ca'）
            k_restraint: 约束力常数 (kJ/(mol·nm²))（默认 10000.0）
            use_com: 是否使用 COM 约束（默认 True）
            platform: 计算平台（CUDA 或 CPU），默认为 config 中的设置

        Returns:
            预平衡后的结构文件路径，如果无 MDP 组件则返回 None
        """
        # 如果未指定 platform，使用 config 中的值（默认为 CUDA）
        if platform is None:
            platform = self.config.simulation.platform

        return self._run_pre_equilibration(
            gpu_id=gpu_id,
            steps=steps,
            mapping=mapping,
            k_restraint=k_restraint,
            use_com=use_com,
            platform=platform,
        )

    def get_pre_equilibrated_structure(self) -> Optional[str]:
        """
        获取预平衡后的结构文件路径

        Returns:
            预平衡结构文件路径，如果未运行预平衡则返回 None
        """
        if self.output_dir is None:
            return None

        preequil_pdb = os.path.join(self.output_dir, 'preequil_final.pdb')
        if os.path.exists(preequil_pdb):
            return preequil_pdb
        return None


# =============================================================================
# Module Exports
# =============================================================================

__all__ = [
    # Configuration
    'CGSimulationConfig',
    'CGComponent',
    'ComponentType',
    'TopologyType',
    'ComputePlatform',
    'SimulationParams',

    # Result
    'SimulationResult',

    # Topology Info
    'TopologyInfo',

    # Simulator
    'CGSimulator',

    # PDB Tools (chain labeling)
    'ChainLabel',
]

