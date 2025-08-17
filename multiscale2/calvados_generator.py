#!/usr/bin/env python3
"""
CALVADOS脚本生成器
根据我们的配置文件自动生成CALVADOS的prepare.py脚本
"""

import os
import yaml
import shutil
from pathlib import Path

class CalvadosGenerator:
    """CALVADOS脚本生成器"""
    
    def __init__(self, config_path):
        """初始化生成器"""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # 保存config_path供其他方法使用
        self.config_path = config_path
        
        # 项目根目录应该是config.yaml的父目录
        self.config_dir = os.path.dirname(os.path.abspath(config_path))
        self.project_root = os.path.dirname(self.config_dir)
        
        # 如果配置文件在子目录中（如example/IDR/），需要调整路径
        if 'example' in config_path and os.path.basename(self.config_dir) in ['IDR', 'MDP']:
            self.project_root = os.path.dirname(self.project_root)
        
        self.calvados_config = self.config.get('cg_calvados', {})
        
        # 读取蛋白质配置
        self.protein_config = self.config.get('protein', {})
        self.protein_name = self.protein_config.get('name', 'protein')
        self.protein_nmol = self.protein_config.get('nmol', 20)
        
        print(f"项目根目录: {self.project_root}")
        print(f"配置文件路径: {config_path}")
        print(f"蛋白质名称: {self.protein_name}")
        print(f"蛋白质数量: {self.protein_nmol}")
        
    def _get_fasta_path(self):
        """获取FASTA文件的正确路径"""
        fasta_relative = self.config['input_files']['sequence_fasta']
        if os.path.isabs(fasta_relative):
            return fasta_relative
        
        # 如果配置文件在子目录中，需要调整路径
        config_dir = os.path.dirname(os.path.abspath(self.config_path))
        if os.path.basename(config_dir) in ['IDR', 'MDP']:
            return os.path.join(config_dir, fasta_relative)
        else:
            return os.path.join(self.project_root, fasta_relative)
        
    def generate_idr_prepare_script(self, output_dir):
        """生成IDR任务的prepare.py脚本"""
        task_type = self.calvados_config.get('task_type', 'IDR')
        
        # 确定residues文件
        if task_type == 'MDP':
            residues_file = 'residues_CALVADOS3.csv'
        else:
            residues_file = 'residues_CALVADOS2.csv'
        
        # 使用配置的蛋白质名称，如果FASTA中存在则验证
        fasta_path = self._get_fasta_path()
        available_names = self._get_available_protein_names_from_fasta(fasta_path)
        
        if self.protein_name not in available_names:
            print(f"⚠️  警告: 配置的蛋白质名称 '{self.protein_name}' 不在FASTA文件中")
            print(f"FASTA文件中的可用序列: {available_names}")
            if available_names:
                self.protein_name = available_names[0]
                print(f"使用第一个可用序列: {self.protein_name}")
            else:
                raise ValueError("FASTA文件中没有找到任何序列")
        else:
            print(f"✓ 找到配置的蛋白质序列: {self.protein_name}")
        
        # 生成prepare.py脚本内容
        script_content = f'''#!/usr/bin/env python3
"""
自动生成的CALVADOS prepare.py脚本
任务类型: {task_type}
"""

import os
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
from argparse import ArgumentParser
from Bio import SeqIO

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',required=True,type=str)
parser.add_argument('--gpu_id',nargs='?',required=True,type=int)
parser.add_argument('--replica',nargs='?',required=True,type=int)
args = parser.parse_args()

cwd = os.getcwd()
N_save = int({self.calvados_config.get('wfreq', 1000)})
N_frames = {self.calvados_config.get('steps', 100000)} // N_save

sysname = f'{{args.name:s}}_{{args.replica:d}}'
residues_file = f'{{cwd}}/input/{residues_file}'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = {self.calvados_config.get('box', [25.0, 25.0, 30.0])}, # nm
  temp = {self.calvados_config.get('temp', 310.0)},
  ionic = {self.calvados_config.get('ionic', 0.15)}, # molar
  pH = {self.calvados_config.get('pH', 7.0)},
  topol = '{self.calvados_config.get('topol', 'random')}',
  slab_width = {self.calvados_config.get('slab_width', 20)},  # slab拓扑的中央高密度区域厚度
  friction = {self.calvados_config.get('friction', 0.01)},

  # INPUT
  ffasta = f'{{cwd}}/input/{os.path.basename(self.config['input_files']['sequence_fasta'])}', # input fasta file
  fresidues = residues_file, # residue definitions

  # RUNTIME SETTINGS
  gpu_id = args.gpu_id,
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  steps = N_frames*N_save, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = '{self.calvados_config.get('platform', 'CUDA')}', # 'CUDA'
  restart = 'checkpoint',
  frestart = '{self.calvados_config.get('frestart', 'restart.chk')}',
  verbose = {str(self.calvados_config.get('verbose', True))},
  slab_eq = {str(self.calvados_config.get('slab_eq', False))},
  steps_eq = 1000,  # 固定为1000步平衡
)

# PATH
path = f'{{cwd}}/{{sysname}}'
output_path = f'{{path}}/data'
subprocess.run(f'mkdir -p {{path}}',shell=True)
subprocess.run(f'mkdir -p {{output_path}}',shell=True)

analyses = f"""
from calvados.analysis import SlabAnalysis, calc_com_traj, calc_contact_map

slab = SlabAnalysis(name="{{sysname:s}}", input_path="{{path:s}}",
                    output_path="{{output_path:s}}", ref_name="{{sysname:s}}", verbose=True)

slab.center(start=400, center_target='all')
slab.calc_profiles()
slab.calc_concentrations()
print(slab.df_results)
slab.plot_density_profiles()
calc_com_traj(path="{{path:s}}",sysname="{{sysname:s}}",output_path="{{output_path:s}}",residues_file="{{residues_file:s}}")
calc_contact_map(path="{{path:s}}",sysname="{{sysname:s}}",output_path="{{output_path:s}}",is_slab=True)
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = {str(self.calvados_config.get('restraint', False))}, # apply restraints
  charge_termini = 'both', # charge N or C or both

  # INPUT
  ffasta = f'{{cwd}}/input/{os.path.basename(self.config['input_files']['sequence_fasta'])}', # input fasta file
  fresidues = f'{{cwd}}/input/{residues_file}', # residue definitions
)

components.add(name=args.name, nmol={self.protein_nmol})  # 使用配置的蛋白质数量

components.write(path,name='components.yaml')
'''
        
        # 写入prepare.py文件
        prepare_path = os.path.join(output_dir, 'prepare.py')
        with open(prepare_path, 'w') as f:
            f.write(script_content)
        
        print(f"✓ 生成IDR prepare.py脚本: {prepare_path}")
        return prepare_path
    
    def generate_mdp_prepare_script(self, output_dir):
        """生成MDP任务的prepare.py脚本"""
        # 确定residues文件
        residues_file = 'residues_CALVADOS3.csv'
        
        # 检查domains文件
        domains_file_path = self.calvados_config.get('fdomains')
        if not domains_file_path:
            raise ValueError("MDP任务需要在cg_calvados配置中指定fdomains")
        
        # 确保domains文件路径是相对于项目根目录的
        if not os.path.isabs(domains_file_path):
            domains_file_path = os.path.join(self.project_root, self.config_dir, domains_file_path)

        if not os.path.exists(domains_file_path):
            raise FileNotFoundError(f"找不到domains文件: {domains_file_path}")

        # 获取PDB文件路径
        pdb_file_name = self.config.get('input_files', {}).get('structure_pdb')
        if not pdb_file_name:
            raise ValueError("MDP任务需要在input_files配置中指定structure_pdb")

        # 生成prepare.py脚本内容
        script_content = f'''#!/usr/bin/env python3
"""
自动生成的CALVADOS prepare.py脚本
任务类型: MDP
"""

import os
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
from argparse import ArgumentParser
from Bio import SeqIO

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',required=True,type=str)
parser.add_argument('--replica',nargs='?',required=True,type=int)
args = parser.parse_args()

cwd = os.getcwd()
N_save = int({self.calvados_config.get('wfreq', 1000)})

sysname = f'{{args.name:s}}_{{args.replica:d}}'
residues_file = f'{{cwd}}/input/{residues_file}'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = {self.calvados_config.get('box', [25.0, 25.0, 30.0])}, # nm
  temp = {self.calvados_config.get('temp', 310.0)},
  ionic = {self.calvados_config.get('ionic', 0.15)}, # molar
  pH = {self.calvados_config.get('pH', 7.0)},
  topol = '{self.calvados_config.get('topol', 'random')}',
  slab_width = {self.calvados_config.get('slab_width', 40)},  # slab拓扑的中央高密度区域厚度
  friction_coeff = {self.calvados_config.get('friction', 0.01)},

  # RUNTIME SETTINGS
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  steps = {self.calvados_config.get('steps', 100000)}, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = '{self.calvados_config.get('platform', 'CUDA')}', # 'CUDA'
  restart = 'checkpoint',
  frestart = '{self.calvados_config.get('frestart', 'restart.chk')}',
  verbose = {str(self.calvados_config.get('verbose', True))},
  slab_eq = {str(self.calvados_config.get('slab_eq', False))},
  steps_eq = 1000,  # 固定为1000步平衡
)

# PATH
path = f'{{cwd}}/{{sysname}}'
subprocess.run(f'mkdir -p {{path}}',shell=True)
subprocess.run(f'mkdir -p data',shell=True)

analyses = f"""
from calvados.analysis import SlabAnalysis, calc_com_traj, calc_contact_map

slab = SlabAnalysis(name="{{sysname:s}}", input_path="{{path:s}}",
  output_path="data", ref_name="{{sysname:s}}", verbose=True)

slab.center(start=0, center_target='all')
slab.calc_profiles()
slab.calc_concentrations()
print(slab.df_results)
slab.plot_density_profiles()

chainid_dict = dict({{sysname:s}} = (0,{self.protein_nmol-1}))
calc_com_traj(path="{{path:s}}",sysname="{{sysname:s}}",output_path="data",residues_file="{{residues_file:s}}",chainid_dict=chainid_dict)
calc_contact_map(path="{{path:s}}",sysname="{{sysname:s}}",output_path="data",chainid_dict=chainid_dict,is_slab=True)
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = {str(self.calvados_config.get('restraint', True))}, # apply restraints
  charge_termini = 'both', # charge N or C or both or none

  # INPUT
  fresidues = residues_file, # residue definitions
  fdomains = f'{{cwd}}/input/{os.path.basename(domains_file_path)}', # domain definitions (harmonic restraints)
  pdb_folder = f'{{cwd}}/input', # directory for pdb and PAE files
)

components.add(
    name='{self.protein_name}', 
    nmol={self.protein_nmol},
    fpdb=f'{{cwd}}/input/{os.path.basename(pdb_file_name)}',
    fdomains=f'{{cwd}}/input/{os.path.basename(domains_file_path)}',
    restraint={str(self.calvados_config.get('restraint', True))}
)
components.write(path,name='components.yaml')
'''
        
        # 写入prepare.py文件
        prepare_path = os.path.join(output_dir, 'prepare.py')
        with open(prepare_path, 'w') as f:
            f.write(script_content)
        
        print(f"✓ 生成MDP prepare.py脚本: {prepare_path}")
        return prepare_path
    
    def _get_protein_name_from_fasta(self, fasta_path):
        """从FASTA文件中获取蛋白质名称"""
        try:
            with open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        return line[1:].strip()
        except:
            pass
        return 'protein'
    
    def _get_available_protein_names_from_fasta(self, fasta_path):
        """从FASTA文件中获取所有可用的蛋白质名称"""
        names = []
        try:
            with open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        name = line[1:].strip()
                        names.append(name)
        except Exception as e:
            print(f"读取FASTA文件时出错: {e}")
        return names
    
    def setup_calvados_environment(self, output_dir):
        """设置CALVADOS运行环境"""
        # 创建input目录
        input_dir = os.path.join(output_dir, 'input')
        os.makedirs(input_dir, exist_ok=True)
        
        # 复制所需文件到输入目录
        if self.calvados_config.get('task_type', 'IDR') == "IDR":
            # 复制FASTA文件
            fasta_path = self._get_fasta_path()
            shutil.copy(fasta_path, os.path.join(input_dir, os.path.basename(fasta_path)))
        elif self.calvados_config.get('task_type', 'IDR') == "MDP":
            # 复制PDB文件
            pdb_file_name = self.config.get('input_files', {}).get('structure_pdb')
            if not pdb_file_name:
                raise ValueError("MDP任务未在config.yaml中指定structure_pdb")
            
            pdb_path = os.path.join(self.project_root, self.config_dir, pdb_file_name)
            if not os.path.exists(pdb_path):
                raise FileNotFoundError(f"找不到PDB文件: {pdb_path}")
            shutil.copy(pdb_path, os.path.join(input_dir, pdb_file_name))

            # 复制domains.yaml文件
            domains_file_path = self.calvados_config.get('fdomains')
            if not domains_file_path:
                raise ValueError("MDP任务未在config.yaml中指定fdomains")

            if not os.path.isabs(domains_file_path):
                domains_file_path = os.path.join(self.project_root, self.config_dir, domains_file_path)
            
            shutil.copy(domains_file_path, os.path.join(input_dir, os.path.basename(domains_file_path)))
        
        # 复制residues.csv文件
        task_type = self.calvados_config.get('task_type', 'IDR')
        residues_suffix = "3" if task_type == "MDP" else "2"
        residues_file_name = f'residues_CALVADOS{residues_suffix}.csv'
        
        residues_src = os.path.join(os.path.dirname(__file__), 'calvados_data', residues_file_name)
        residues_dst = os.path.join(input_dir, residues_file_name)
        shutil.copy(residues_src, residues_dst)
        
        print(f"✓ 设置CALVADOS环境: {input_dir}")
        return input_dir
    
    def generate_and_run(self, output_dir, protein_name=None, gpu_id=0, replica=1):
        """生成脚本并运行CALVADOS"""
        # 如果没有指定蛋白质名称，使用配置中的名称
        if protein_name is None:
            protein_name = self.protein_name
        
        # 设置环境
        self.setup_calvados_environment(output_dir)
        
        # 生成prepare.py脚本
        task_type = self.calvados_config.get('task_type', 'IDR')
        if task_type == 'MDP':
            prepare_path = self.generate_mdp_prepare_script(output_dir)
        else:
            prepare_path = self.generate_idr_prepare_script(output_dir)
        
        # 运行prepare.py
        print(f"\n运行CALVADOS prepare.py...")
        print(f"蛋白质: {protein_name}")
        print(f"蛋白质数量: {self.protein_nmol}")
        print(f"GPU ID: {gpu_id}")
        print(f"副本: {replica}")
        
        # 构建命令
        if task_type == 'MDP':
            cmd = f"cd {output_dir} && python prepare.py --name {protein_name} --replica {replica}"
        else:
            cmd = f"cd {output_dir} && python prepare.py --name {protein_name} --gpu_id {gpu_id} --replica {replica}"
        
        print(f"执行命令: {cmd}")
        return cmd
