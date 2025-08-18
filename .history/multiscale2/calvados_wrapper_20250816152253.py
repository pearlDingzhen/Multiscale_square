import os
import shutil
import yaml
import logging
from typing import Dict, Any
from .calvados import sim, cfg

logger = logging.getLogger(__name__)

class CalvadosWrapperError(Exception):
    """CALVADOS包装器错误异常"""
    pass

def run_calvados_simulation(config: Dict[str, Any], output_dir: str, components_path: str):
    """
    A wrapper to run the CALVADOS simulation as the first stage.
    It prepares a compatible config file from the main config.
    
    Args:
        config: Main configuration dictionary
        output_dir: Directory to run the simulation in
        components_path: Path to the components.yaml file
        
    Raises:
        CalvadosWrapperError: If simulation setup or execution fails
    """
    logger.info("Setting up CALVADOS simulation...")
    
    try:
        # 验证输入参数
        if not config:
            raise CalvadosWrapperError("Configuration dictionary is empty")
        
        if not os.path.exists(output_dir):
            raise CalvadosWrapperError(f"Output directory does not exist: {output_dir}")
        
        if not os.path.exists(components_path):
            raise CalvadosWrapperError(f"Components file does not exist: {components_path}")
        
        # 提取CALVADOS配置
        calvados_config_dict = config.get('cg_calvados', {})
        if not calvados_config_dict:
            raise CalvadosWrapperError("Configuration for 'cg_calvados' not found in config file.")
        
        logger.info("Extracted CALVADOS configuration parameters")
        
        # 验证CALVADOS配置的完整性
        _validate_calvados_config_dict(calvados_config_dict)
        
        # 创建临时配置文件
        temp_config_path = os.path.join(output_dir, 'cg_config.yaml')
        _create_calvados_config_file(calvados_config_dict, temp_config_path)
        
        # 复制components文件
        temp_components_path = os.path.join(output_dir, 'components.yaml')
        _copy_components_file(components_path, temp_components_path)
        
        # 运行CALVADOS模拟
        logger.info("Starting CALVADOS simulation...")
        sim.run(path=output_dir, fconfig='cg_config.yaml', fcomponents='components.yaml')
        
        logger.info("CALVADOS simulation completed successfully")
        
    except Exception as e:
        logger.error(f"CALVADOS simulation failed: {e}")
        raise CalvadosWrapperError(f"CALVADOS simulation failed: {e}")

def _validate_calvados_config_dict(config_dict: Dict[str, Any]):
    """验证CALVADOS配置字典的完整性"""
    logger.info("Validating CALVADOS configuration dictionary...")
    
    # 检查必需的参数
    required_params = [
        'platform', 'temp', 'box', 'steps', 'sysname', 
        'eps_lj', 'ionic', 'pH', 'wfreq', 'logfreq'
    ]
    
    missing_params = [param for param in required_params if param not in config_dict]
    if missing_params:
        raise CalvadosWrapperError(f"Missing required CALVADOS parameters: {missing_params}")
    
    # 验证参数值
    if config_dict['temp'] <= 0:
        raise CalvadosWrapperError("Temperature must be positive")
    
    if len(config_dict['box']) != 3:
        raise CalvadosWrapperError("Box must be a list of 3 dimensions")
    
    if config_dict['steps'] <= 0:
        raise CalvadosWrapperError("Number of steps must be positive")
    
    if config_dict['eps_lj'] <= 0:
        raise CalvadosWrapperError("Lennard-Jones epsilon must be positive")
    
    if not (0 <= config_dict['pH'] <= 14):
        raise CalvadosWrapperError("pH must be between 0 and 14")
    
    logger.info("CALVADOS configuration validation passed")

def _create_calvados_config_file(config_dict: Dict[str, Any], config_path: str):
    """创建CALVADOS配置文件"""
    try:
        with open(config_path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, indent=2)
        
        logger.info(f"Created CALVADOS config file: {config_path}")
        
    except Exception as e:
        raise CalvadosWrapperError(f"Failed to create CALVADOS config file: {e}")

def _copy_components_file(source_path: str, dest_path: str):
    """复制components文件"""
    try:
        shutil.copy2(source_path, dest_path)
        logger.info(f"Copied components file: {source_path} -> {dest_path}")
        
    except Exception as e:
        raise CalvadosWrapperError(f"Failed to copy components file: {e}")