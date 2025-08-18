#!/usr/bin/env python3
"""
CALVADOS Script Generator
Automatically generates CALVADOS prepare.py scripts based on our configuration files
"""

import os
import yaml
import shutil
from pathlib import Path

class CalvadosGenerator:
    """CALVADOS Script Generator"""
    
    def __init__(self, config_path):
        """Initialize the generator"""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # Save config_path for use by other methods
        self.config_path = config_path
        self.config_dir = os.path.dirname(os.path.abspath(config_path))
        
        self.calvados_config = self.config.get('cg_calvados', {})
        
        # Prefer top-level task_type, fallback to cg_calvados.task_type, default to IDP
        self.task_type = (
            self.config.get('task_type')
            or self.calvados_config.get('task_type')
            or 'IDP'
        )
        
        # Read protein configuration
        self.protein_config = self.config.get('protein', {})
        self.protein_name = self.protein_config.get('name', 'protein')
        self.protein_nmol = self.protein_config.get('nmol', 20)
        
        print(f"Configuration file path: {config_path}")
        print(f"Protein name: {self.protein_name}")
        print(f"Number of protein molecules: {self.protein_nmol}")
        
    def _get_fasta_path(self):
        """Get the correct path to the FASTA file"""
        fasta_relative = self.config['input_files']['sequence_fasta']
        if os.path.isabs(fasta_relative):
            return fasta_relative
        
        return os.path.join(self.config_dir, fasta_relative)
        
    def generate_idp_prepare_script(self, output_dir, gpu_id=0):
        """Generate prepare.py script for IDP tasks"""
        task_type = self.task_type
        
        # Determine residues file
        if task_type == 'MDP':
            residues_file = 'residues_CALVADOS3.csv'
        else:
            residues_file = 'residues_CALVADOS2.csv'
        
        # Use configured protein name, verify if it exists in FASTA
        fasta_path = self._get_fasta_path()
        available_names = self._get_available_protein_names_from_fasta(fasta_path)
        
        if self.protein_name not in available_names:
            print(f"⚠️  Warning: Configured protein name '{self.protein_name}' not found in FASTA file")
            print(f"Available sequences in FASTA file: {available_names}")
            if available_names:
                self.protein_name = available_names[0]
                print(f"Using first available sequence: {self.protein_name}")
            else:
                raise ValueError("No sequences found in FASTA file")
        else:
            print(f"✓ Found configured protein sequence: {self.protein_name}")
        
        # Generate prepare.py script content
        script_content = f'''#!/usr/bin/env python3
"""
Auto-generated CALVADOS prepare.py script
Task type: {task_type}
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
  slab_width = {self.calvados_config.get('slab_width', 20)},  # thickness of central high-density region for slab topology
  friction = {self.calvados_config.get('friction', 0.01)},

  # INPUT
  ffasta = f'{{cwd}}/input/{os.path.basename(self.config['input_files']['sequence_fasta'])}', # input fasta file
  fresidues = residues_file, # residue definitions

  # RUNTIME SETTINGS
  gpu_id = {gpu_id},  # Use the gpu_id parameter passed to generate_and_run
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  steps = N_frames*N_save, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = '{self.calvados_config.get('platform', 'CUDA')}', # 'CUDA'
  restart = 'checkpoint',
  frestart = '{self.calvados_config.get('frestart', 'restart.chk')}',
  verbose = {str(self.calvados_config.get('verbose', True))},
  slab_eq = {str(self.calvados_config.get('slab_eq', False))},
  steps_eq = 1000,  # fixed to 1000 steps for equilibration
)

# PATH
path = f'{{cwd}}/{{sysname}}'
output_path = f'{{path}}/data'
subprocess.run(f'mkdir -p {{path}}',shell=True)
subprocess.run(f'mkdir -p {{output_path}}',shell=True)

analyses = ""

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

components.add(name=args.name, nmol={self.protein_nmol})  # use configured number of protein molecules

components.write(path,name='components.yaml')
'''
        
        # Write prepare.py file
        prepare_path = os.path.join(output_dir, 'prepare.py')
        with open(prepare_path, 'w') as f:
            f.write(script_content)
        
        print(f"✓ Generated IDP prepare.py script: {prepare_path}")
        return prepare_path
    
    def generate_mdp_prepare_script(self, output_dir):
        """Generate prepare.py script for MDP tasks"""
        # Determine residues file
        residues_file = 'residues_CALVADOS3.csv'
        
        # Check domains file
        domains_file_path = self.calvados_config.get('fdomains')
        if not domains_file_path:
            raise ValueError("MDP tasks require fdomains to be specified in cg_calvados configuration")
        
        # Ensure domains file path is relative to project root
        if not os.path.isabs(domains_file_path):
            domains_file_path = os.path.join(self.config_dir, domains_file_path)

        if not os.path.exists(domains_file_path):
            raise FileNotFoundError(f"Domains file not found: {domains_file_path}")

        # Get PDB file path
        pdb_file_name = self.config.get('input_files', {}).get('structure_pdb')
        if not pdb_file_name:
            raise ValueError("MDP tasks require structure_pdb to be specified in input_files configuration")

        # Generate prepare.py script content
        script_content = f'''#!/usr/bin/env python3
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
  slab_width = {self.calvados_config.get('slab_width', 40)},  # thickness of central high-density region for slab topology
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
  steps_eq = 1000,  # fixed to 1000 steps for equilibration
)

# PATH
path = f'{{cwd}}/{{sysname}}'
subprocess.run(f'mkdir -p {{path}}',shell=True)
subprocess.run(f'mkdir -p data',shell=True)

analyses = ""

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
        
        # Write prepare.py file
        prepare_path = os.path.join(output_dir, 'prepare.py')
        with open(prepare_path, 'w') as f:
            f.write(script_content)
        
        print(f"✓ Generated MDP prepare.py script: {prepare_path}")
        return prepare_path
    
    def _get_protein_name_from_fasta(self, fasta_path):
        """Get protein name from FASTA file"""
        try:
            with open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        return line[1:].strip()
        except:
            pass
        return 'protein'
    
    def _get_available_protein_names_from_fasta(self, fasta_path):
        """Get all available protein names from FASTA file"""
        names = []
        try:
            with open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        name = line[1:].strip()
                        names.append(name)
        except Exception as e:
            print(f"Error reading FASTA file: {e}")
        return names
    
    def setup_calvados_environment(self, output_dir):
        """Set up CALVADOS simulation environment"""
        # Create input directory
        input_dir = os.path.join(output_dir, 'input')
        os.makedirs(input_dir, exist_ok=True)
        
        # Copy required files to input directory
        if self.task_type == "IDP":
            # Copy FASTA file
            fasta_path = self._get_fasta_path()
            shutil.copy(fasta_path, os.path.join(input_dir, os.path.basename(fasta_path)))
        elif self.task_type == "MDP":
            # Copy PDB file
            pdb_file_name = self.config.get('input_files', {}).get('structure_pdb')
            if not pdb_file_name:
                raise ValueError("structure_pdb not specified in config.yaml for MDP task")
            
            pdb_path = os.path.join(self.config_dir, pdb_file_name)
            if not os.path.exists(pdb_path):
                raise FileNotFoundError(f"PDB file not found: {pdb_path}")
            shutil.copy(pdb_path, os.path.join(input_dir, pdb_file_name))

            # Copy domains.yaml file
            domains_file_path = self.calvados_config.get('fdomains')
            if not domains_file_path:
                raise ValueError("fdomains not specified in config.yaml for MDP task")

            if not os.path.isabs(domains_file_path):
                domains_file_path = os.path.join(self.config_dir, domains_file_path)
            
            shutil.copy(domains_file_path, os.path.join(input_dir, os.path.basename(domains_file_path)))
        
        # Copy residues.csv file (use global task_type)
        task_type = self.task_type
        residues_suffix = "3" if task_type == "MDP" else "2"
        residues_file_name = f'residues_CALVADOS{residues_suffix}.csv'
        
        residues_src = os.path.join(os.path.dirname(__file__), 'calvados_data', residues_file_name)
        residues_dst = os.path.join(input_dir, residues_file_name)
        shutil.copy(residues_src, residues_dst)
        
        print(f"✓ Set up CALVADOS environment: {input_dir}")
        return input_dir
    
    def generate_and_run(self, output_dir, protein_name=None, gpu_id=0, replica=1):
        """Generate scripts and run CALVADOS"""
        # If protein name is not specified, use the name from configuration
        if protein_name is None:
            protein_name = self.protein_name
        
        # Set up environment
        self.setup_calvados_environment(output_dir)
        
        # Generate prepare.py script
        task_type = self.task_type
        if task_type == 'MDP':
            prepare_path = self.generate_mdp_prepare_script(output_dir)
        else:
            prepare_path = self.generate_idp_prepare_script(output_dir, gpu_id)
        
        # Run prepare.py
        print(f"\nRunning CALVADOS prepare.py...")
        print(f"Protein: {protein_name}")
        print(f"Number of protein molecules: {self.protein_nmol}")
        print(f"GPU ID: {gpu_id}")
        print(f"Replica: {replica}")
        
        # Build prepare command (Step 1)
        prepare_cmd = f"cd {output_dir} && python prepare.py --name {protein_name} --replica {replica}"
        
        print(f"Step 1 - Prepare command: {prepare_cmd}")
        
        # Build run command (Step 2)
        sysname = f"{protein_name}_{replica}"
        run_cmd = f"cd {output_dir} && python {sysname}/run.py --path {sysname}"
        
        print(f"Step 2 - Run command: {run_cmd}")
        
        # Return both commands as a list
        return [prepare_cmd, run_cmd]
