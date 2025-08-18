
import os
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
from argparse import ArgumentParser
from Bio import SeqIO

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',required=True,type=str)
parser.add_argument('--replica',nargs='?',required=True,type=int)
parser.add_argument('--temp', nargs='?',required=True,type=float)
parser.add_argument('--number', nargs='?',required=True,type=int)
parser.add_argument('--box', nargs='?',required=True,type=int)
args = parser.parse_args()

cwd = os.getcwd()
N_save = int(5e4) # 50000 step = 500 ps = 1 ns 

sysname = f'{args.name:s}_{args.replica:d}'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [args.box, args.box, args.box], # nm
  temp = args.temp, # K
  ionic = 0.1, # molar
  pH = 7.5,
  topol = 'grid',
  friction_coeff = 0.001,
  

  # RUNTIME SETTINGS
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  steps = 500*N_save, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CUDA', # 'CUDA'
  restart = 'checkpoint',
  frestart = 'restart.chk',
  verbose = True
)

# PATH
path = f'{cwd}/{sysname}'
subprocess.run(f'mkdir -p {path}',shell=True)
subprocess.run(f'mkdir -p data',shell=True)


with open(f'{path}/system_info.txt','w') as f:
    f.write(f'Replica: {args.replica:d}\n')
    f.write(f'Number of molecules: {args.number:d}\n')
    f.write(f'Box size: {args.box:d} nm\n')
    f.write(f'Temperature: {args.temp:.2f} K\n')



analyses = f"""
from calvados.analysis import calc_slab_profiles

calc_slab_profiles(path="{path:s}",name="{sysname:s}",output_folder="data",ref_atoms="all",start=0)
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = args.number, # number of molecules
  restraint = True, # apply restraints
  charge_termini = 'both', # charge N or C or both or none

  # INPUT
  fresidues = f'{cwd}/input/residues_CALVADOS3.csv', # residue definitions
  fdomains = f'{cwd}/input/domains.yaml', # domain definitions (harmonic restraints)
  pdb_folder = f'{cwd}/input', # directory for pdb and PAE files

)

components.add(name=args.name, nmol=args.number)
components.write(path,name='components.yaml')
