#!/usr/bin/env python3
"""
Backmap Module

Handles conversion from coarse-grained (CG) structures to all-atom (AA) representations.
Supports two input modes:
1. ms2 cg output: Reads final.pdb from CG simulation output
2. user provided: Standardizes PDB using calvados and then backmaps
"""

import os
import tempfile
from pathlib import Path
from typing import Optional, List, Dict
from enum import Enum
from dataclasses import dataclass

from openmm import Vec3
import openmm as mm
import openmm.unit as unit
from openmm.app import PDBFile

from .cg import CGSimulationConfig, ComponentType, BackmapConfig
from .calvados_wrapper import CalvadosWrapper


class SourceType(Enum):
    """输入源类型"""
    MS2_CG = "ms2_cg"
    USER_PROVIDED = "user_provided"


@dataclass
class PreparedInput:
    """准备好的输入数据"""
    pdb_path: str
    model_type: str  # "ResidueBasedModel" 或 "CalphaBasedModel"
    config: Optional[CGSimulationConfig] = None  # 原始配置（如果有）


@dataclass
class BackmapResult:
    """Backmap 结果"""
    success: bool
    output_pdb: str
    input_pdb: str
    model_type: str
    errors: List[str] = None
    
    def __post_init__(self):
        if self.errors is None:
            self.errors = []


def standardize_pdb_with_calvados(pdb_path: str, config: CGSimulationConfig, output_pdb: str) -> str:
    """
    使用 calvados 构建 system 并标准化 PDB 格式（CA + 三字母代码）
    
    这个方法从 cg.py 的 _convert_mpipi_pdb_to_calvados_format 提取而来，
    可以用于任何需要标准化的 PDB 文件。
    
    Args:
        pdb_path: 输入的 PDB 文件路径
        config: CGSimulationConfig 对象（包含 components 信息）
        output_pdb: 输出 PDB 文件路径
        
    Returns:
        输出 PDB 文件路径
    """
    temp_dir = tempfile.mkdtemp(prefix='calvados_topology_')
    try:
        # 1. 使用 calvados wrapper 构建 system（这会生成 top.pdb）
        wrapper = CalvadosWrapper(config)
        
        # 写入配置文件到临时目录
        wrapper._write_to_dir(temp_dir, gpu_id=0, verbose=False)
        
        # 创建 calvados Sim 对象并构建 system
        from multiscale2.extern.ms2_calvados.calvados import sim as calvados_sim
        from yaml import safe_load
        
        with open(f'{temp_dir}/config.yaml', 'r') as stream:
            calvados_config = safe_load(stream)
        with open(f'{temp_dir}/components.yaml', 'r') as stream:
            components = safe_load(stream)
        
        # 创建 Sim 对象（不运行模拟，只构建 system）
        calvados_sim_obj = calvados_sim.Sim(temp_dir, calvados_config, components)
        calvados_sim_obj.build_system()
        
        # 2. 从 calvados 的 top.pdb 读取 topology
        calvados_top_pdb = os.path.join(temp_dir, 'top.pdb')
        if not os.path.exists(calvados_top_pdb):
            raise FileNotFoundError(f"Calvados top.pdb not found: {calvados_top_pdb}")
        
        calvados_pdb_file = PDBFile(calvados_top_pdb)
        calvados_topology = calvados_pdb_file.topology
        
        # 3. 从输入 PDB 读取坐标
        input_pdb_file = PDBFile(pdb_path)
        input_positions = input_pdb_file.positions
        
        # 验证原子数匹配
        n_calvados = calvados_topology.getNumAtoms()
        n_input = len(input_positions)
        
        if n_calvados != n_input:
            raise ValueError(
                f"原子数不匹配: calvados={n_calvados}, input={n_input}. "
                f"请检查 component 定义是否一致。"
            )
        
        # 4. 转换坐标单位：Angstrom -> nm（PDB 文件使用 Angstrom，OpenMM 使用 nm）
        positions_nm_quantity = input_positions.in_units_of(unit.nanometer)
        
        # 5. 获取盒子向量（从输入 PDB 的 topology 或使用 config）
        box_vectors = None
        try:
            box_vectors = input_pdb_file.topology.getPeriodicBoxVectors()
            if box_vectors is None:
                raise ValueError("No box vectors in input topology")
        except:
            # 如果获取失败，使用 config 的 box 创建盒子向量
            pass
        
        # 如果没有从 PDB 获取到，使用 config 的 box 创建盒子向量
        if box_vectors is None:
            box = config.box
            box_vectors = [
                mm.Vec3(box[0], 0, 0) * unit.nanometer,
                mm.Vec3(0, box[1], 0) * unit.nanometer,
                mm.Vec3(0, 0, box[2]) * unit.nanometer
            ]
        
        # 设置盒子向量到 topology
        calvados_topology.setPeriodicBoxVectors(box_vectors)
        
        # 6. 使用 calvados topology + 输入坐标保存 PDB
        with open(output_pdb, 'w') as f:
            PDBFile.writeFile(
                calvados_topology,
                positions_nm_quantity,
                f,
                keepIds=False  # 不保留原始 ID，让 OpenMM 重新编号
            )
        
        return output_pdb
        
    finally:
        # 清理临时目录
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)


class BackmapSimulator:
    """Backmap 模拟器，处理 CG 到 AA 的转换"""
    
    def __init__(self, config: Optional[CGSimulationConfig] = None, backmap_config: Optional[BackmapConfig] = None):
        """
        初始化 backmap 配置
        
        Args:
            config: CGSimulationConfig（可选，用于 user provided 模式）
            backmap_config: BackmapConfig（可选，从 config.backmap 或 CLI 参数获取）
        """
        self.config = config
        self.backmap_config = backmap_config or BackmapConfig()
        
        # 如果 config 中有 backmap 配置，合并它（CLI 参数优先）
        if config and config.backmap:
            if not backmap_config or not backmap_config.model_type:
                if config.backmap.model_type:
                    self.backmap_config.model_type = config.backmap.model_type
            if not backmap_config or not backmap_config.device:
                if config.backmap.device:
                    self.backmap_config.device = config.backmap.device
            if not backmap_config or not backmap_config.output_dir:
                if config.backmap.output_dir:
                    self.backmap_config.output_dir = config.backmap.output_dir
    
    def detect_source_type(self, input_path: str) -> SourceType:
        """
        检测输入是 ms2 cg 输出还是 user provided
        
        Args:
            input_path: 输入路径（目录或文件）
            
        Returns:
            SourceType.MS2_CG 或 SourceType.USER_PROVIDED
        """
        path = Path(input_path)
        
        # 检查是否是目录（ms2 cg 输出）
        if path.is_dir():
            # 检查目录结构特征
            has_final_pdb = (path / "final.pdb").exists()
            has_simulation_log = (path / "simulation.log").exists()
            
            if has_final_pdb or has_simulation_log:
                return SourceType.MS2_CG
        
        # 检查是否是 PDB 文件（user provided）
        if path.is_file() and path.suffix.lower() == '.pdb':
            return SourceType.USER_PROVIDED
        
        raise ValueError(f"Cannot determine source type for: {input_path}")
    
    def find_config_yaml(self, cg_output_dir: str, system_name: str, explicit_config: Optional[str] = None) -> Optional[str]:
        """
        查找 config.yaml 文件
        
        Args:
            cg_output_dir: CG 输出目录
            system_name: 系统名称
            explicit_config: 用户显式提供的 config 路径
            
        Returns:
            config.yaml 路径，如果找不到返回 None
        """
        # 优先使用用户显式提供的
        if explicit_config:
            config_path = Path(explicit_config)
            if config_path.exists():
                return str(config_path.resolve())
            else:
                raise FileNotFoundError(f"Config file not found: {explicit_config}")
        
        # 备选：查找当前工作目录
        pwd_config = Path.cwd() / f"{system_name}.yaml"
        if pwd_config.exists():
            return str(pwd_config)
        
        # 备选：查找父目录
        parent_config = Path(cg_output_dir).parent / f"{system_name}.yaml"
        if parent_config.exists():
            return str(parent_config)
        
        return None
    
    def select_model_type(self, config: CGSimulationConfig, force_field: Optional[str] = None) -> str:
        """
        选择 CG model 类型
        
        规则：
        - 如果用户指定了 model_type，使用用户指定的
        - 否则：Calvados + MDP → ResidueBasedModel，其他 → CalalphaBasedModel
        
        Args:
            config: CGSimulationConfig
            force_field: 力场名称（可选，用于检测）
            
        Returns:
            "ResidueBasedModel" 或 "CalphaBasedModel"
        """
        # 如果用户指定了 model_type，使用用户指定的
        if self.backmap_config.model_type and self.backmap_config.model_type != 'auto':
            return self.backmap_config.model_type
        
        # 自动选择逻辑
        has_mdp = any(c.type == ComponentType.MDP for c in config.components)
        
        # 检测力场（如果未提供）
        if force_field is None:
            force_field = self._detect_force_field_from_dir(config)
        
        is_calvados = force_field == 'calvados'
        
        if is_calvados and has_mdp:
            return "ResidueBasedModel"
        else:
            return "CalphaBasedModel"
    
    def _detect_force_field_from_dir(self, config: CGSimulationConfig) -> Optional[str]:
        """
        从配置或目录结构检测力场
        
        这是一个简单的启发式方法，可以通过检查输出目录结构来推断力场
        """
        # 如果 config 有相关信息，可以在这里添加检测逻辑
        # 目前返回 None，让调用者提供
        return None
    
    def prepare_ms2_cg_input(self, cg_output_dir: str, config_path: Optional[str] = None) -> PreparedInput:
        """
        准备 ms2 cg 输出用于 backmap
        
        Args:
            cg_output_dir: CG 输出目录路径
            config_path: 显式提供的 config.yaml 路径（可选）
            
        Returns:
            PreparedInput 对象
        """
        cg_output_path = Path(cg_output_dir)
        
        # 1. 读取 final.pdb
        final_pdb = cg_output_path / "final.pdb"
        if not final_pdb.exists():
            raise FileNotFoundError(f"final.pdb not found in {cg_output_dir}")
        
        # 2. 查找并加载 config.yaml
        # 从目录名推断 system_name（假设格式为 {system_name}_CG）
        system_name = cg_output_path.name.replace('_CG', '')
        
        config_yaml_path = self.find_config_yaml(str(cg_output_dir), system_name, config_path)
        
        if config_yaml_path:
            config = CGSimulationConfig.from_yaml(config_yaml_path)
            # 更新 simulator 的 config（用于后续使用）
            self.config = config
            # 合并 backmap 配置（CLI 参数优先）
            if config.backmap:
                if not self.backmap_config.model_type and config.backmap.model_type:
                    self.backmap_config.model_type = config.backmap.model_type
                if not self.backmap_config.device and config.backmap.device:
                    self.backmap_config.device = config.backmap.device
                if not self.backmap_config.output_dir and config.backmap.output_dir:
                    self.backmap_config.output_dir = config.backmap.output_dir
        else:
            # 如果没有找到 config，创建一个最小配置（仅用于 model type 选择）
            # 这种情况下，我们无法准确判断是否有 MDP，默认使用 CalalphaBasedModel
            config = None
        
        # 3. 检测力场（通过目录结构）
        force_field = self._detect_force_field_from_directory(cg_output_path)
        
        # 4. 选择 CG model
        if config:
            model_type = self.select_model_type(config, force_field)
        else:
            # 如果没有 config，默认使用 CalalphaBasedModel
            model_type = "CalphaBasedModel"
        
        return PreparedInput(
            pdb_path=str(final_pdb),
            model_type=model_type,
            config=config
        )
    
    def _detect_force_field_from_directory(self, cg_output_path: Path) -> Optional[str]:
        """从目录结构检测力场类型"""
        # 检查子目录名称
        subdirs = [d.name for d in cg_output_path.iterdir() if d.is_dir()]
        
        if 'Mpipi-Recharged' in subdirs or 'mpipi_recharged' in subdirs:
            return 'mpipi_recharged'
        elif 'HPS' in subdirs or 'hps' in subdirs:
            return 'hps_urry'
        elif 'COCOMO' in subdirs or 'cocomo' in subdirs:
            return 'cocomo'
        elif 'raw' in subdirs:
            # raw 目录通常表示 calvados
            return 'calvados'
        
        return None
    
    def prepare_user_provided_input(self, pdb_path: str, config: CGSimulationConfig) -> PreparedInput:
        """
        准备 user provided PDB 用于 backmap
        
        Args:
            pdb_path: 用户提供的 PDB 文件路径
            config: CGSimulationConfig（必须包含 components 信息）
            
        Returns:
            PreparedInput 对象
        """
        if not config or not config.components:
            raise ValueError("User provided mode requires config with components")
        
        # 1. 使用 calvados 构建标准化 PDB
        temp_standardized = tempfile.mktemp(suffix='.pdb', prefix='standardized_')
        standardized_pdb = standardize_pdb_with_calvados(pdb_path, config, temp_standardized)
        
        # 2. 选择 model type
        model_type = self.select_model_type(config, force_field='calvados')
        
        return PreparedInput(
            pdb_path=standardized_pdb,
            model_type=model_type,
            config=config
        )
    
    def run(self, input_path: str, config_path: Optional[str] = None, output_dir: Optional[str] = None) -> BackmapResult:
        """
        执行 backmap
        
        Args:
            input_path: 输入路径（CG 输出目录或 PDB 文件）
            config_path: 配置文件路径（可选）
            output_dir: 输出目录（可选，默认 {system_name}_backmap）
            
        Returns:
            BackmapResult 对象
        """
        result = BackmapResult(
            success=False,
            output_pdb="",
            input_pdb=input_path,
            model_type="",
            errors=[]
        )
        
        try:
            # 1. 检测输入类型
            source_type = self.detect_source_type(input_path)
            
            # 2. 准备输入
            if source_type == SourceType.MS2_CG:
                prepared = self.prepare_ms2_cg_input(input_path, config_path)
            else:  # USER_PROVIDED
                if not config_path:
                    raise ValueError("User provided mode requires config.yaml via -f option")
                config = CGSimulationConfig.from_yaml(config_path)
                # 更新 simulator 的 config
                self.config = config
                # 合并 backmap 配置（CLI 参数优先）
                if config.backmap:
                    if not self.backmap_config.model_type and config.backmap.model_type:
                        self.backmap_config.model_type = config.backmap.model_type
                    if not self.backmap_config.device and config.backmap.device:
                        self.backmap_config.device = config.backmap.device
                    if not self.backmap_config.output_dir and config.backmap.output_dir:
                        self.backmap_config.output_dir = config.backmap.output_dir
                prepared = self.prepare_user_provided_input(input_path, config)
            
            result.model_type = prepared.model_type
            
            # 3. 确定输出目录（优先级：CLI 参数 > config.backmap.output_dir > 默认）
            if output_dir is None:
                if self.backmap_config.output_dir:
                    output_dir = self.backmap_config.output_dir
                elif prepared.config:
                    output_dir = f"{prepared.config.system_name}_backmap"
                else:
                    # 从输入路径推断
                    input_name = Path(input_path).stem
                    if input_name.endswith('_CG'):
                        input_name = input_name[:-3]
                    output_dir = f"{input_name}_backmap"
            
            os.makedirs(output_dir, exist_ok=True)
            
            # 4. 确定输出文件名
            output_pdb = os.path.join(output_dir, f"{Path(prepared.pdb_path).stem}.aa.pdb")
            
            # 5. 执行 backmap
            from multiscale2.extern.ms2_cg2all import convert_cg2all
            
            device = self.backmap_config.device
            convert_cg2all(
                in_pdb_fn=prepared.pdb_path,
                out_fn=output_pdb,
                model_type=prepared.model_type,
                fix_atom=False,
                device=device,
                write_ssbond=False  # 默认不写入二硫键记录
            )
            
            result.success = True
            result.output_pdb = output_pdb
            
        except Exception as e:
            result.errors.append(str(e))
            result.success = False
        
        return result
