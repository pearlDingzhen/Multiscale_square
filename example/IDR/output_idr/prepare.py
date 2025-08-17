#!/usr/bin/env python3
"""
自动生成的CALVADOS prepare.py脚本
任务类型: IDR
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
N_save = int(10000)
N_frames = 10000000 // N_save

sysname = f'{args.name:s}_{args.replica:d}'
residues_file = f'{cwd}/input/residues_CALVADOS2.csv'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [10.0, 10.0, 40.0], # nm
  temp = 310.0,
  ionic = 0.15, # molar
  pH = 7.0,
  topol = 'slab',
  slab_width = 20,  # slab拓扑的中央高密度区域厚度
  friction = 0.01,

  # INPUT
  ffasta = f'{cwd}/input/protein.fasta', # input fasta file
  fresidues = residues_file, # residue definitions

  # RUNTIME SETTINGS
  gpu_id = args.gpu_id,
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  steps = N_frames*N_save, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CUDA', # 'CUDA'
  restart = 'checkpoint',
  frestart = 'restart.chk',
  verbose = True,
  slab_eq = False,
  steps_eq = 1000,  # 固定为1000步平衡
)

# PATH
path = f'{cwd}/{sysname}'
output_path = f'{path}/data'
subprocess.run(f'mkdir -p {path}',shell=True)
subprocess.run(f'mkdir -p {output_path}',shell=True)

analyses = f"""
from calvados.analysis import SlabAnalysis, calc_com_traj, calc_contact_map

slab = SlabAnalysis(name="{sysname:s}", input_path="{path:s}",
                    output_path="{output_path:s}", ref_name="{sysname:s}", verbose=True)

slab.center(start=400, center_target='all')
slab.calc_profiles()
slab.calc_concentrations()
print(slab.df_results)
slab.plot_density_profiles()
calc_com_traj(path="{path:s}",sysname="{sysname:s}",output_path="{output_path:s}",residues_file="{residues_file:s}")
calc_contact_map(path="{path:s}",sysname="{sysname:s}",output_path="{output_path:s}",is_slab=True)
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = False, # apply restraints
  charge_termini = 'both', # charge N or C or both

  # INPUT
  ffasta = f'{cwd}/input/protein.fasta', # input fasta file
  fresidues = f'{cwd}/input/residues_CALVADOS2.csv', # residue definitions
)

components.add(name=args.name, nmol=20)  # 使用配置的蛋白质数量

components.write(path,name='components.yaml')
