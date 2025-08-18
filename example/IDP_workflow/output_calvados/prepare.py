#!/usr/bin/env python3
"""
Auto-generated CALVADOS prepare.py script
Task type: IDP
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
N_frames = 200000 // N_save

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
  slab_width = 20,  # thickness of central high-density region for slab topology
  friction = 0.01,

  # INPUT
  ffasta = f'{cwd}/input/protein.fasta', # input fasta file
  fresidues = residues_file, # residue definitions

  # RUNTIME SETTINGS
  gpu_id = 0,  # Use the gpu_id parameter passed to generate_and_run
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  steps = N_frames*N_save, # number of simulation steps
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
output_path = f'{path}/data'
subprocess.run(f'mkdir -p {path}',shell=True)
subprocess.run(f'mkdir -p {output_path}',shell=True)

analyses = ""

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

components.add(name=args.name, nmol=20)  # use configured number of protein molecules

components.write(path,name='components.yaml')
