import os
import shutil
import yaml

# 尝试导入系统安装的CALVADOS
try:
    import calvados
    from calvados import sim
    print("✓ 使用系统安装的CALVADOS")
except ImportError:
    # 如果系统没有安装，回退到本地版本
    try:
        from .calvados import sim
        print("⚠️  使用本地CALVADOS版本")
    except ImportError:
        raise ImportError("CALVADOS未安装，请安装CALVADOS: pip install calvados")

def run_calvados_simulation(config, output_dir, components_path):
    """
    A wrapper to run the CALVADOS simulation as the first stage.
    It prepares a compatible config file from the main config.
    """
    # Extract the 'cg_calvados' section for compatibility
    calvados_config_dict = config.get('cg_calvados', {})
    if not calvados_config_dict:
        raise ValueError("Configuration for 'cg_calvados' not found in config file.")

    # 确定任务类型
    task_type = calvados_config_dict.get('task_type', 'IDR')  # 默认为IDR
    print(f"✓ 任务类型: {task_type}")
    
    # 根据任务类型选择residues文件
    if task_type == 'MDP':
        residues_file = 'residues_CALVADOS3.csv'
        print("✓ 使用MDP residues文件 (CALVADOS3)")
    else:  # IDR
        residues_file = 'residues_CALVADOS2.csv'
        print("✓ 使用IDR residues文件 (CALVADOS2)")
    
    # 获取residues文件的绝对路径
    current_dir = os.path.dirname(os.path.abspath(__file__))
    residues_path = os.path.join(current_dir, 'calvados_data', residues_file)
    
    if not os.path.exists(residues_path):
        raise FileNotFoundError(f"Residues文件未找到: {residues_path}")

    # 添加CALVADOS必需的默认参数
    default_params = {
        'slab_eq': False,
        'bilayer_eq': False,
        'ext_force': False,
        'restart': None,
        'k_eq': 1000.0,
        'ext_force_expr': '',
        'friction': 0.01,
        'verbose': True,
    }
    
    # 根据任务类型添加特定参数
    if task_type == 'MDP':
        default_params.update({
            'restraint': True,
            'fdomains': None,  # 需要用户提供
        })
    else:  # IDR
        default_params.update({
            'restraint': False,
        })
    
    # 合并配置，确保所有必需参数都存在
    for key, default_value in default_params.items():
        if key not in calvados_config_dict:
            calvados_config_dict[key] = default_value

    # Write a temporary, CALVADOS-compatible config file
    temp_config_path = os.path.join(output_dir, 'cg_config.yaml')
    with open(temp_config_path, 'w') as f:
        yaml.dump(calvados_config_dict, f)

    # 读取并修改components文件
    with open(components_path, 'r') as f:
        components_data = yaml.safe_load(f)
    
    # 更新residues文件路径
    if 'defaults' in components_data:
        components_data['defaults']['fresidues'] = residues_path
        print(f"✓ 更新residues文件路径: {residues_path}")
    
    # 如果是MDP任务，需要添加domains文件
    if task_type == 'MDP':
        domains_file = calvados_config_dict.get('fdomains')
        if domains_file and os.path.exists(domains_file):
            # 复制domains文件到输出目录
            temp_domains_path = os.path.join(output_dir, 'domains.yaml')
            shutil.copy(domains_file, temp_domains_path)
            components_data['defaults']['fdomains'] = 'domains.yaml'
            print(f"✓ 复制domains文件: {temp_domains_path}")
        else:
            print("⚠️  MDP任务需要domains文件，但未找到或未指定")
    
    # 写入修改后的components文件
    temp_components_path = os.path.join(output_dir, 'components.yaml')
    with open(temp_components_path, 'w') as f:
        yaml.dump(components_data, f)

    print(f"✓ CALVADOS配置文件已创建: {temp_config_path}")
    print(f"✓ Components文件已创建: {temp_components_path}")
    
    # 显示配置内容
    print("\nCALVADOS配置内容:")
    with open(temp_config_path, 'r') as f:
        print(f.read())

    # Run the simulation using the CALVADOS entry point
    print("\n开始运行CALVADOS模拟...")
    sim.run(path=output_dir, fconfig='cg_config.yaml', fcomponents='components.yaml')