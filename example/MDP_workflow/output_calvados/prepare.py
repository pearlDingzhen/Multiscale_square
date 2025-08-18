#!/usr/bin/env python3
"""
Auto-generated CALVADOS prepare.py script
Task type: MDP
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
N_save = int(5000)

sysname = f'{args.name:s}_{args.replica:d}'
residues_file = f'{cwd}/input/residues_CALVADOS3.csv'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [10.0, 10.0, 40.0], # nm
  temp = 310.0,
  ionic = 0.15, # molar
  pH = 7.5,
  topol = 'slab',
  slab_width = 40,  # thickness of central high-density region for slab topology
  friction_coeff = 0.01,

  # RUNTIME SETTINGS
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  steps = 500000, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CUDA', # 'CUDA'
  restart = 'checkpoint',
  frestart = 'restart.chk',
  verbose = True,
  slab_eq = False,
  steps_eq = 1000,  # fixed to 1000 steps for equilibration
)

# PATH
path = f'{cwd}/{sysname}'
subprocess.run(f'mkdir -p {path}',shell=True)
subprocess.run(f'mkdir -p data',shell=True)

analyses = ""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = True, # apply restraints
  charge_termini = 'both', # charge N or C or both or none

  # INPUT
  fresidues = residues_file, # residue definitions
  fdomains = f'{cwd}/input/domains.yaml', # domain definitions (harmonic restraints)
  pdb_folder = f'{cwd}/input', # directory for pdb and PAE files
)

components.add(
    name='TDP43', 
    nmol=15,
    fpdb=f'{cwd}/input/TDP43.pdb',
    fdomains=f'{cwd}/input/domains.yaml',
    restraint=True
)
components.write(path,name='components.yaml')
