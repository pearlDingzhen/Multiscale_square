#!/usr/bin/env python3
"""
Stage 3: OpenMM Optimization
Auto-generated script for generating topology and optimizing structure with OpenMM.
"""

import os
import sys
import yaml
import subprocess
import shutil
import glob
from pathlib import Path
import openmm
# ============================================================================
# USER CONFIGURABLE VARIABLES
# ============================================================================
# These variables can be modified by the user if needed

# Input directory from Stage 2 (leave None to auto-detect)
INPUT_DIR = None

# Input structure file (leave "LATEST" to auto-detect)
STRUCTURE_FILE = "LATEST"

# Output directory for OpenMM optimization
OUTPUT_DIR = "output_openmm"

# Topology generation settings
TOPOLOGY_SETTINGS = {
    "force_regenerate": False,  # Set to True to regenerate topology even if exists
}

# OpenMM optimization settings
OPENMM_SETTINGS = {
    "temperature": 325,         # K
    "pressure": 1,             # bar
    "time_step": 0.004,        # ps
    "platform": "CUDA",        # CUDA, OpenCL, CPU
    "precision": "double"      # single, mixed, double
}

# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

def load_config():
    """Load configuration from YAML file."""
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def auto_detect_input_dir():
    """Auto-detect input directory from Stage 2 output."""
    backmap_dir = "output_backmap"
    if os.path.exists(backmap_dir):
        return backmap_dir
    
    print(f"Default backmap directory not found: {backmap_dir}")
    return None

def read_fasta_sequence(fasta_path):
    """Read sequence from FASTA file."""
    sequences = {}
    if not os.path.exists(fasta_path):
        return sequences
    
    name = None
    seq_chunks = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = ''.join(seq_chunks)
                name = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if name is not None:
        sequences[name] = ''.join(seq_chunks)
    return sequences

def extract_sequence_from_pdb(pdb_path):
    """Extract sequence from PDB file."""
    import MDAnalysis as mda
    u = mda.Universe(pdb_path)
    
    # Get protein residues and extract sequence
    protein = u.select_atoms("protein")
    residues = protein.residues
    
    # Map 3-letter to 1-letter codes
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    sequence = ""
    for residue in residues:
        resname = residue.resname.strip()
        if resname in three_to_one:
            sequence += three_to_one[resname]
        else:
            print(f"Warning: Unknown residue {resname}, using 'X'")
            sequence += 'X'
    
    return sequence

def get_protein_sequence():
    """Get protein sequence from configuration."""
    config = load_config()
    task_type = config.get('task_type', 'IDP')
    
    if task_type == 'IDP':
        # Read from FASTA file
        input_files = config.get('input_files', {})
        fasta_file = input_files.get('sequence_fasta')
        if not fasta_file or not os.path.exists(fasta_file):
            print(f"Error: FASTA file not found: {fasta_file}")
            sys.exit(1)
        
        sequences = read_fasta_sequence(fasta_file)
        protein_name = config.get('protein', {}).get('name')
        
        if protein_name and protein_name in sequences:
            return sequences[protein_name]
        elif len(sequences) == 1:
            return list(sequences.values())[0]
        else:
            print(f"Error: Could not determine sequence from FASTA file")
            sys.exit(1)
            
    elif task_type == 'MDP':
        # Read from PDB file
        input_files = config.get('input_files', {})
        pdb_file = input_files.get('structure_pdb')
        if not pdb_file or not os.path.exists(pdb_file):
            print(f"Error: PDB file not found: {pdb_file}")
            sys.exit(1)
        
        return extract_sequence_from_pdb(pdb_file)
    
    else:
        print(f"Error: Unknown task type: {task_type}")
        sys.exit(1)

def generate_topology(sequence, output_dir, protein_name=None):
    """Generate PACE topology using prepare_peptide.py."""
    print("="*60)
    print("Step 1: Generating PACE Topology")
    print("="*60)
    
    print(f"Sequence: {sequence}")
    print(f"Sequence length: {len(sequence)}")
    if protein_name:
        print(f"Protein name: {protein_name}")
    
    # Use protein name if provided, otherwise use first 10 characters of sequence
    if protein_name:
        safe_name = protein_name
    else:
        safe_name = sequence[:10] if len(sequence) > 10 else sequence
    
    # Path to prepare_peptide.py using multiscale2 package path
    import multiscale2
    pace_builder_dir = Path(multiscale2.__path__[0]) / "pace_top_builder"
    if not pace_builder_dir.exists():
        print(f"Error: PACE top builder directory not found: {pace_builder_dir}")
        sys.exit(1)
    
    prepare_script = pace_builder_dir / "prepare_peptide.py"
    if not prepare_script.exists():
        print(f"Error: prepare_peptide.py not found: {prepare_script}")
        sys.exit(1)
    
    # Change to PACE builder directory and run prepare_peptide.py
    original_cwd = os.getcwd()
    # IMPORTANT: Resolve topology_dir to an absolute path BEFORE changing directory
    topology_dir = Path(output_dir).resolve() / "topology"
    topology_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        os.chdir(pace_builder_dir)
        
        # Run prepare_peptide.py without --output-dir, use default out/ directory
        cmd = [sys.executable, str(prepare_script), sequence]
        if protein_name:
            cmd.extend(["--name", safe_name])
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error running prepare_peptide.py:")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            sys.exit(1)
        
        print("✓ Topology generation completed successfully!")
        print(result.stdout)
        
        # Move output files from pace_builder_dir/out to the absolute topology_dir and cleanup
        out_dir = pace_builder_dir / "out"  # This is an absolute path
        print(f"Looking for output files in: {out_dir}")
        
        if out_dir.exists():
            # Define source and destination files
            source_top = out_dir / "topol.top"
            dest_top = topology_dir / "PACE.top" # Rename to PACE.top
            source_pdb = out_dir / "system.pdb"
            dest_pdb = topology_dir / "system.pdb"

            # Move topology file
            if source_top.exists():
                shutil.move(str(source_top), str(dest_top))
                print(f"Moved: {source_top} -> {dest_top}")
            else:
                print(f"Warning: topol.top not found in {out_dir}")

            # Move structure file
            if source_pdb.exists():
                shutil.move(str(source_pdb), str(dest_pdb))
                print(f"Moved: {source_pdb} -> {dest_pdb}")
            else:
                print(f"Warning: system.pdb not found in {out_dir}")
            
            # Remove the now-empty out directory
            shutil.rmtree(out_dir)
            print(f"Cleaned up: Removed {out_dir} directory")
        else:
            print(f"Warning: Output directory not found after running prepare_peptide.py: {out_dir}")
        
        return dest_pdb, dest_top # Return absolute paths to the final files
        
    finally:
        os.chdir(original_cwd)

def modify_topology_nmol(topol_path, nmol):
    """Modify topology file to match number of molecules from CG simulation."""
    print(f"\nModifying topology for {nmol} molecules...")
    
    with open(topol_path, 'r') as f:
        lines = f.readlines()
    
    # Find and modify the molecules section
    in_molecules = False
    for i, line in enumerate(lines):
        if line.strip().startswith('[ molecules ]'):
            in_molecules = True
            continue
        if in_molecules and line.strip() and not line.startswith(';'):
            parts = line.strip().split()
            if len(parts) >= 2:
                # Modify the number of molecules
                parts[1] = str(nmol)
                # Ensure a single newline character at the end
                lines[i] = ' '.join(parts) + '\n'
                break
    
    # Write back to file
    with open(topol_path, 'w') as f:
        f.writelines(lines)
    
    print(f"✓ Updated topology to include {nmol} molecules")

def run_pdb2gmx(input_structure, output_dir):
    """Run pdb2gmx on the backmapped structure."""
    print("="*60)  
    print("Step 2: Processing Structure with pdb2gmx")
    print("="*60)
    
    # Import gromacs wrapper
    try:
        import gromacs
    except ImportError:
        print("Error: GromacsWrapper is not installed")
        print("Please install it using: pip install GromacsWrapper")
        sys.exit(1)
    
    # Define absolute paths. All operations will be based on the script's launch directory.
    cwd = Path.cwd()
    gromacs_dir = (cwd / output_dir / "gromacs").resolve()
    gromacs_dir.mkdir(parents=True, exist_ok=True)
    
    input_pdb_path = (cwd / input_structure).resolve()
    if not input_pdb_path.exists():
        print(f"FATAL: Input PDB file not found at resolved absolute path: {input_pdb_path}")
        sys.exit(1)

    # Define output file paths
    local_pdb_path = gromacs_dir / input_pdb_path.name
    output_gro_path = gromacs_dir / "processed.gro"
    output_top_path = gromacs_dir / "topol.top"

    try:
        # Copy input structure to the gromacs directory
        shutil.copy2(str(input_pdb_path), str(local_pdb_path))
        print(f"Copied input structure to: {local_pdb_path}")
        
        print(f"Processing structure: {local_pdb_path}")
        print("Using PACE force field...")
        
        # Run pdb2gmx using absolute paths, without changing directory
        gromacs.pdb2gmx(
            f=str(local_pdb_path),
            ff="pace-new", 
            water="no",
            ignh=True,
            o=str(output_gro_path),  # Generate .gro format for OpenMM
            p=str(output_top_path),
            i="posre.itp"  # Generate position restraints
        )
        
        print("✓ pdb2gmx completed successfully!")
        print(f"  - Output structure: {output_gro_path}")
        print(f"  - Output topology (unused): {output_top_path}")
        
        return output_gro_path, output_top_path
        
    except Exception as e:
        print(f"Error running pdb2gmx: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def add_ring_constraints_to_system(system, topology, 
                                    dihedral_k=-50.0, 
                                    angle_k=400.0,
                                    ring_dihedral_dict=None,
                                    ring_angle_dict=None,
                                    ring_diagonal_dict=None):
	"""
	为芳香族氨基酸的环结构添加完整的约束：二面角、键角和对角线距离限制
	
	Parameters
	----------
	system : openmm.System
		OpenMM系统对象
	topology : openmm.app.Topology
		拓扑对象
	dihedral_k : float, default=-50.0
		二面角力常数 (kJ/mol)
	angle_k : float, default=400.0
		键角力常数 (kJ/mol/rad^2)
	ring_dihedral_dict : dict, optional
		环二面角字典，如果为None则使用默认值
	ring_angle_dict : dict, optional
		环键角字典，如果为None则使用默认值
	ring_diagonal_dict : dict, optional
		环对角线距离限制字典，如果为None则使用默认值
		
	Returns
	-------
	forces : dict
		包含添加的所有力的字典
	"""
	from openmm import PeriodicTorsionForce, HarmonicAngleForce, CustomBondForce
	import openmm.unit as unit
	import math
	
	degToRad = math.pi / 180.0
	
	# 默认二面角字典
	if ring_dihedral_dict is None:
		ring_dihedral_dict = {
			'HIS': {'CD2-CG-ND1-CE1': 0,
			        'CE1-NE2-CD2-CG': 0,
			        'CG-ND1-CE1-NE2': 0,
			        'ND1-CE1-NE2-CD2': 0,
			        'NE2-CD2-CG-ND1': 0},
			'PHE': {'CD1-CE1-CZ-CE2': 0,
			        'CD2-CG-CD1-CE1': 0,
			        'CE1-CZ-CE2-CD2': 0,
			        'CE2-CD2-CG-CD1': 0,
			        'CG-CD1-CE1-CZ': 0,
			        'CZ-CE2-CD2-CG': 0},
			'TRP': {'CD1-CG-CD2-CE3': 180,
			        'CD1-NE1-CE2-CD2': 0,
			        'CD1-NE1-CE2-CZ2': 180,
			        'CD2-CE2-CZ2-CH2': 0,
			        'CD2-CG-CD1-NE1': 0,
			        'CE2-CD2-CG-CD1': 0,
			        'CE2-CZ2-CH2-CZ3': 0,
			        'CE3-CD2-CE2-CZ2': 0,
			        'CG-CD1-NE1-CE2': 0,
			        'CG-CD2-CE3-CZ3': 180,
			        'CH2-CZ3-CE3-CD2': 0,
			        'CZ2-CE2-CD2-CG': 180,
			        'CZ2-CH2-CZ3-CE3': 0,
			        'CZ3-CE3-CD2-CE2': 0,
			        'NE1-CE2-CD2-CE3': 180,
			        'NE1-CE2-CD2-CG': 0,
			        'NE1-CE2-CZ2-CH2': 180},
			'TYR': {'CD1-CE1-CZ-CE2': 0,
			        'CD2-CG-CD1-CE1': 0,
			        'CE1-CZ-CE2-CD2': 0,
			        'CE2-CD2-CG-CD1': 0,
			        'CG-CD1-CE1-CZ': 0,
			        'CZ-CE2-CD2-CG': 0}
		}
	
	# 默认键角字典
	if ring_angle_dict is None:
		ring_angle_dict = {
			'HIS': {'ring': {
				'ND1-CE1-NE2': 108.0,
				'CE1-NE2-CD2': 108.0,
				'NE2-CD2-CG': 108.0,
				'CD2-CG-ND1': 108.0,
				'CG-ND1-CE1': 108.0,
				'ND1-CD2-CG': 108.0,
				'CD2-ND1-CE1': 108.0,
				'CD2-CE1-NE2': 108.0
			}},
			'PHE': {'ring': {
				'CG-CD1-CE1': 120.0,
				'CD1-CE1-CZ': 120.0,
				'CE1-CZ-CE2': 120.0,
				'CZ-CE2-CD2': 120.0,
				'CE2-CD2-CG': 120.0,
				'CD2-CG-CD1': 120.0,
				'CD1-CD2-CE2': 120.0,
				'CD1-CZ-CE2': 120.0,
				'CE1-CD1-CD2': 120.0,
				'CE1-CZ-CD2': 120.0
			}},
			'TYR': {'ring': {
				'CG-CD1-CE1': 120.0,
				'CD1-CE1-CZ': 120.0,
				'CE1-CZ-CE2': 120.0,
				'CZ-CE2-CD2': 120.0,
				'CE2-CD2-CG': 120.0,
				'CD2-CG-CD1': 120.0,
				'CD1-CD2-CE2': 120.0,
				'CD1-CZ-CE2': 120.0,
				'CE1-CD1-CD2': 120.0,
				'CE1-CZ-CD2': 120.0
			}},
			'TRP': {
				'five_membered': {
					'CG-CD1-NE1': 108.0,
					'CD1-NE1-CE2': 108.0,
					'NE1-CE2-CD2': 108.0,
					'CE2-CD2-CG': 108.0,
					'CD2-CG-CD1': 108.0,
					'CD1-CE2-NE1': 108.0,
					'CD1-CD2-CE2': 108.0,
					'NE1-CD2-CG': 108.0
				},
				'six_membered': {
					'CG-CD2-CE3': 120.0,
					'CD2-CE3-CZ3': 120.0,
					'CE3-CZ3-CH2': 120.0,
					'CZ3-CH2-CZ2': 120.0,
					'CH2-CZ2-CE2': 120.0,
					'CZ2-CE2-CD2': 120.0,
					'CD2-CZ3-CE3': 120.0,
					'CD2-CH2-CZ2': 120.0,
					'CD2-CE3-CZ2': 120.0,
					'CE3-CD2-CZ2': 120.0,
					'CE3-CH2-CZ3': 120.0
				}
			}
		}
	
	# 默认对角线距离限制字典 (Type 10 Distance Restraints)
	if ring_diagonal_dict is None:
		ring_diagonal_dict = {
			'PHE': {
				'ring': {
					('CG', 'CZ'): (0.270, 0.290, 0.310, 5000),
					('CD1', 'CE2'): (0.270, 0.290, 0.310, 5000),
					('CD2', 'CE1'): (0.270, 0.290, 0.310, 5000)
				}
			},
			'TYR': {
				'ring': {
					('CG', 'CZ'): (0.270, 0.290, 0.310, 5000),
					('CD1', 'CE2'): (0.270, 0.290, 0.310, 5000),
					('CD2', 'CE1'): (0.270, 0.290, 0.310, 5000)
				}
			},
			'HIS': {
				'ring': {
					('CG', 'CE1'): (0.215, 0.235, 0.255, 5000),
					('CG', 'NE2'): (0.215, 0.235, 0.255, 5000),
					('ND1', 'NE2'): (0.215, 0.235, 0.255, 5000),
					('ND1', 'CD2'): (0.215, 0.235, 0.255, 5000),
					('CE1', 'CD2'): (0.215, 0.235, 0.255, 5000)
				}
			},
			'TRP': {
				'six_membered': {
					('CD2', 'CH2'): (0.270, 0.290, 0.310, 5000),
					('CE2', 'CZ3'): (0.270, 0.290, 0.310, 5000),
					('CZ2', 'CE3'): (0.270, 0.290, 0.310, 5000)
				},
				'five_membered': {
					('CG', 'NE1'): (0.220, 0.240, 0.260, 5000),
					('CG', 'CE2'): (0.220, 0.240, 0.260, 5000),
					('CD1', 'CE2'): (0.220, 0.240, 0.260, 5000),
					('CD1', 'CD2'): (0.220, 0.240, 0.260, 5000),
					('NE1', 'CD2'): (0.220, 0.240, 0.260, 5000)
				}
			}
		}
	
	# 创建力对象
	dihedral_force = PeriodicTorsionForce()
	angle_force = HarmonicAngleForce()
	
	# Type 10 Distance Restraints 势能函数
	distance_restraint_expression = (
		"step(r0 - r) * 0.5 * kdr * (r - r0)^2 + "
		"step(r - r1) * step(r2 - r) * 0.5 * kdr * (r - r1)^2 + "
		"step(r - r2) * 0.5 * kdr * (r2 - r1) * (2*r - r2 - r1);"
	)
	distance_force = CustomBondForce(distance_restraint_expression)
	distance_force.addPerBondParameter("r0")  # low
	distance_force.addPerBondParameter("r1")  # up1
	distance_force.addPerBondParameter("r2")  # up2
	distance_force.addPerBondParameter("kdr")  # force constant
	
	# 辅助函数：查找原子索引
	def find_atom_index(residue, atom_name):
		for atom in residue.atoms():
			if atom.name == atom_name:
				return atom.index
		return None
	
	# 添加二面角约束
	dihedrals_added = 0
	for residue in topology.residues():
		residue_name = residue.name
		if residue_name in ring_dihedral_dict:
			dihedral_info = ring_dihedral_dict[residue_name]
			for dihedral_name, target_angle in dihedral_info.items():
				atom_names = dihedral_name.split('-')
				if len(atom_names) == 4:
					indices = [find_atom_index(residue, name) for name in atom_names]
					if all(i is not None for i in indices):
						dihedral_force.addTorsion(
							indices[0], indices[1], indices[2], indices[3],
							1,  # periodicity
							target_angle * degToRad,  # phase
							dihedral_k  # k
						)
						dihedrals_added += 1
	
	# 添加键角约束
	angles_added = 0
	for residue in topology.residues():
		residue_name = residue.name
		if residue_name in ring_angle_dict:
			angle_data = ring_angle_dict[residue_name]
			for ring_type, angles in angle_data.items():
				for angle_name, ideal_angle in angles.items():
					atom_names = angle_name.split('-')
					if len(atom_names) == 3:
						indices = [find_atom_index(residue, name) for name in atom_names]
						if all(i is not None for i in indices):
							angle_force.addAngle(
								indices[0], indices[1], indices[2],
								ideal_angle * degToRad,  # theta0
								angle_k  # k
							)
							angles_added += 1
	
	# 添加对角线距离限制 (Type 10 Distance Restraints)
	distances_added = 0
	for residue in topology.residues():
		residue_name = residue.name
		if residue_name in ring_diagonal_dict:
			diagonal_data = ring_diagonal_dict[residue_name]
			for ring_type, diagonal_pairs in diagonal_data.items():
				for (atom1_name, atom2_name), (low, up1, up2, kdr) in diagonal_pairs.items():
					idx1 = find_atom_index(residue, atom1_name)
					idx2 = find_atom_index(residue, atom2_name)
					if idx1 is not None and idx2 is not None:
						distance_force.addBond(
							idx1, idx2,
							[low, up1, up2, kdr]  # r0, r1, r2, kdr
						)
						distances_added += 1
	
	# 将力添加到系统
	forces = {}
	if dihedrals_added > 0:
		system.addForce(dihedral_force)
		forces['dihedral'] = dihedral_force
	if angles_added > 0:
		system.addForce(angle_force)
		forces['angle'] = angle_force
	if distances_added > 0:
		system.addForce(distance_force)
		forces['distance'] = distance_force
	
	return forces


def run_openmm_optimization(structure_gro, topology_top, output_dir, task_type):
    """Run OpenMM optimization with integrated logic."""
    # Ensure all necessary modules are imported within the function scope
    from pathlib import Path
    import os
    import sys
    import shutil
    
    print("="*60)
    print("Step 3: OpenMM Structure Optimization") 
    print("="*60)
    
    # Change to OpenMM directory and run optimization
    openmm_dir = Path(structure_gro).parent  # Files should already be in openmm directory
    original_cwd = os.getcwd()
    try:
        os.chdir(openmm_dir)
        print(f"Working in OpenMM directory: {openmm_dir}")
        
        print("Running OpenMM minimization...")
        
        # Import OpenMM and PACE modules
        import openmm.unit as unit
        import openmm
        from openmm.app import GromacsGroFile, Simulation, PDBFile
        from openmm import LangevinIntegrator, CustomExternalForce, Platform
        from multiscale2 import PACE2openmm as PACE
        
        # Configuration
        conf_name = 'conf.gro'
        top_name = 'PACE.top'
        time_step = 0.004
        ref_t = 325
        ref_p = 1
        epsilon_r = 15
        
        # Platform setup
        properties = {'Precision': 'double'}
        platform = Platform.getPlatformByName("CUDA")
        
        # Load structure and topology
        conf = GromacsGroFile(conf_name)
        gro_position = conf.getPositions()
        box_vectors = conf.getPeriodicBoxVectors()
        
        # Find GROMACS topology directory using environment variables
        import os
        from pathlib import Path
        
        def find_gromacs_topology_dir():
            """
            Finds the GROMACS topology directory by checking standard 
            GROMACS environment variables.
            
            Returns:
                str: The absolute path to the 'top' directory if found, otherwise None.
            """
            # GMXSHARE is the standard for modern GROMACS versions (>= 5.x)
            gmxshare_path = os.environ.get('GMXSHARE')
            if gmxshare_path:
                # The 'top' directory is typically under 'gromacs/' within GMXSHARE
                top_path = Path(gmxshare_path) / 'gromacs' / 'top'
                if top_path.is_dir():
                    print(f"Found GROMACS topology path via GMXSHARE: {top_path}")
                    return str(top_path)
            
            # GMXDATA was used for older versions (4.x)
            gmxdata_path = os.environ.get('GMXDATA')
            if gmxdata_path:
                # GMXDATA often pointed to a directory containing 'top'
                top_path = Path(gmxdata_path) / 'top'
                if top_path.is_dir():
                    print(f"Found GROMACS topology path via GMXDATA: {top_path}")
                    return str(top_path)
            
            print("Warning: Could not find GROMACS topology directory via environment variables.")
            print("Please ensure you have sourced the GMXRC file before running the script.")
            return None
        
        GMXLIB = find_gromacs_topology_dir()
        
        if GMXLIB is None:
            print("Warning: Could not find GROMACS topology directory automatically.")
            print("Please set GMXSHARE or GMXDATA environment variable or ensure GROMACS is properly installed.")
            # Fallback to common locations as last resort
            fallback_candidates = [
                '/usr/share/gromacs/top/',
                '/usr/local/share/gromacs/top/',
                '/opt/gromacs/share/gromacs/top/',
                os.path.expanduser('~/gromacs/share/gromacs/top/'),
                '/mnt/hdd1/tianxj_out/gmx2022.6.cuda/share/gromacs/top/'
            ]
            
            for candidate in fallback_candidates:
                if os.path.exists(candidate):
                    GMXLIB = candidate
                    print(f"Using fallback GROMACS topology path: {GMXLIB}")
                    break
            
            if GMXLIB is None:
                print("Error: No GROMACS topology directory found. Please check your GROMACS installation.")
                sys.exit(1)
        
        top = PACE.PACETopFile(
            top_name,
            periodicBoxVectors=box_vectors,
            defines={},
            epsilon_r=epsilon_r,
            includeDir=GMXLIB
        )
        
        ################################################################################
        # Step 1: Minimize with Gaussian "Bulldozer" Force (Initial untangling)
        ################################################################################
        
        print("Step 1: Minimizing with Gaussian 'Bulldozer' Force (Initial untangling)...")
        
        # Gaussian mode uses pure geometric repulsion to untangle knots
        system_gaussian = top.create_system(
            nonbonded_cutoff=1.1 * unit.nanometer, 
            add_nonbonded_force=True,
            nonbonded_type="gaussian"
        )
        
        # 为高斯模式系统也添加环二面角约束
        print("为高斯模式系统添加芳香族氨基酸环二面角约束...")
        ring_forces_gaussian = add_ring_constraints_to_system(system_gaussian, top.topology, dihedral_k=-50.0, angle_k=400.0)
        
        # Add CA restraint only for IDP tasks
        if task_type == 'IDP':
            restraint_gaussian = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
            system_gaussian.addForce(restraint_gaussian)
            restraint_gaussian.addGlobalParameter('k', 100*unit.kilojoules_per_mole/unit.nanometer)
            restraint_gaussian.addPerParticleParameter('x0')
            restraint_gaussian.addPerParticleParameter('y0')
            restraint_gaussian.addPerParticleParameter('z0')
            
            for atom in top.topology.atoms():
                if atom.name == 'CA':
                    restraint_gaussian.addParticle(atom.index, gro_position[atom.index])
        
        integrator_gaussian = LangevinIntegrator(ref_t * unit.kelvin,
                                               10.0 / unit.picosecond,
                                               time_step * unit.picosecond)
        integrator_gaussian.setRandomNumberSeed(0)
        
        simulation_gaussian = Simulation(top.topology, system_gaussian, integrator_gaussian, platform, properties)
        simulation_gaussian.context.setPositions(conf.getPositions())
        
        # Energy before Gaussian minimization
        state = simulation_gaussian.context.getState(getEnergy=True, getForces=False)
        print(f'Energy (Gaussian) before minimization: {state.getPotentialEnergy()}')
        
        # Minimize with Gaussian force
        simulation_gaussian.minimizeEnergy(maxIterations=5000, tolerance=100)
        
        # Save Gaussian minimized structure
        state_gaussian = simulation_gaussian.context.getState(getPositions=True, enforcePeriodicBox=True)
        PDBFile.writeFile(simulation_gaussian.topology, state_gaussian.getPositions(), open('opti1.pdb', 'w'))
        
        state = simulation_gaussian.context.getState(getEnergy=True, getForces=False)
        print(f'Energy (Gaussian) after minimization: {state.getPotentialEnergy()}')
        print('Gaussian minimized structure saved as opti1.pdb')
        
        ################################################################################
        # Step 2: Minimize with Standard Martini Force (Final optimization)
        # TEST MODE: Skipping softcore steps, going directly from Gaussian to Standard
        ################################################################################
        
        print("\nStep 2: Minimizing with Standard Martini Force (Final optimization)...")
        print("TEST MODE: Skipping softcore steps, going directly from Gaussian to Standard")
        
        # Standard mode with full LJ and electrostatics
        system_standard = top.create_system(
            nonbonded_cutoff=1.1 * unit.nanometer,
            add_nonbonded_force=True,
            nonbonded_type="standard"
        )
        
        # 为标准模式系统也添加环二面角约束
        print("为标准模式系统添加芳香族氨基酸环二面角约束...")
        ring_forces_standard = add_ring_constraints_to_system(system_standard, top.topology, dihedral_k=-50.0, angle_k=400.0)
        
        # Add CA restraint only for IDP tasks
        if task_type == 'IDP':
            restraint_standard = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
            system_standard.addForce(restraint_standard)
            restraint_standard.addGlobalParameter('k', 100*unit.kilojoules_per_mole/unit.nanometer)
            restraint_standard.addPerParticleParameter('x0')
            restraint_standard.addPerParticleParameter('y0')
            restraint_standard.addPerParticleParameter('z0')
            
            for atom in top.topology.atoms():
                if atom.name == 'CA':
                    restraint_standard.addParticle(atom.index, gro_position[atom.index])
        
        integrator_standard = LangevinIntegrator(ref_t * unit.kelvin,
                                                10.0 / unit.picosecond,
                                                time_step * unit.picosecond)
        integrator_standard.setRandomNumberSeed(0)
        
        simulation_standard = Simulation(top.topology, system_standard, integrator_standard, platform, properties)
        # Use Gaussian minimized positions as starting point
        simulation_standard.context.setPositions(state_gaussian.getPositions())
        
        # Energy before standard minimization
        state = simulation_standard.context.getState(getEnergy=True, getForces=False)
        print(f'Energy (Standard) before minimization: {state.getPotentialEnergy()}')
        
        # Minimize with standard force
        simulation_standard.minimizeEnergy(maxIterations=5000, tolerance=100)
        
        # Final state
        state_standard = simulation_standard.context.getState(getPositions=True, enforcePeriodicBox=True)
        PDBFile.writeFile(simulation_standard.topology, state_standard.getPositions(), open('opti2.pdb', 'w'))
        
        state = simulation_standard.context.getState(getEnergy=True, getForces=False)
        print(f'Energy (Standard) after minimization: {state.getPotentialEnergy()}')
        
        print('\033[32mMinimization successful! Final structure saved as opti2.pdb\033[0m')
        
        # Check for output files
        output_files = ["opti1.pdb", "opti2.pdb"]
        found_files = []
        for output_file in output_files:
            if os.path.exists(output_file):
                found_files.append(output_file)
                
        if found_files:
            print(f"Generated optimized structures: {found_files}")
            return [openmm_dir / f for f in found_files]
        else:
            print("Warning: No optimized structure files found")
            return []
            
    except Exception as e:
        print(f"Error during OpenMM optimization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        # Clean up .itp files in the current directory before changing back
        try:
            import glob
            itp_files = glob.glob("*.itp")
            if itp_files:
                for itp_file in itp_files:
                    os.remove(itp_file)
                print(f"Cleaned up {len(itp_files)} .itp file(s)")
        except Exception as cleanup_error:
            print(f"Warning: Failed to clean up .itp files: {cleanup_error}")
        
        os.chdir(original_cwd)

def run_openmm_workflow():
    """Execute the complete OpenMM workflow."""
    print("="*60)
    print("Stage 3: OpenMM Optimization Workflow")
    print("="*60)
    
    # Load configuration
    config = load_config()
    
    # Determine input directory
    input_dir = INPUT_DIR or auto_detect_input_dir()
    if not input_dir:
        print("Error: Could not determine input directory from Stage 2")
        print("Please set INPUT_DIR variable or ensure Stage 2 completed successfully")
        sys.exit(1)
    
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Find input structure file
    if STRUCTURE_FILE == "LATEST":
        # Look for final backmapped structure
        candidates = [
            "backmapped_aa_final_ter.pdb",
            "backmapped_aa_final.pdb", 
            "backmapped_aa_box.pdb"
        ]
        
        input_structure = None
        for candidate in candidates:
            candidate_path = os.path.join(input_dir, candidate)
            if os.path.exists(candidate_path):
                input_structure = candidate_path
                break
        
        if not input_structure:
            print(f"Error: No suitable structure file found in {input_dir}")
            print(f"Looked for: {candidates}")
            sys.exit(1)
    else:
        input_structure = os.path.join(input_dir, STRUCTURE_FILE)
        if not os.path.exists(input_structure):
            print(f"Error: Structure file not found: {input_structure}")
            sys.exit(1)
    
    print(f"Using structure: {os.path.basename(input_structure)}")
    
    # Step 1: Generate topology
    sequence = get_protein_sequence()
    config = load_config()
    protein_name = config.get('protein', {}).get('name', None)
    topology_pdb, topology_top = generate_topology(sequence, OUTPUT_DIR, protein_name)
    
    # Modify topology for correct number of molecules
    protein_config = config.get('protein', {})
    nmol = protein_config.get('nmol', 1)
    modify_topology_nmol(topology_top, nmol)
    
    # Step 2: Process structure with pdb2gmx to get conf.gro
    processed_pdb, _ = run_pdb2gmx(input_structure, OUTPUT_DIR)
    
    # Step 3: Run OpenMM optimization
    # Copy PACE.top to openmm directory and use it with conf.gro
    openmm_dir = Path(OUTPUT_DIR) / "openmm"
    openmm_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy PACE topology to openmm directory
    pace_top_dest = openmm_dir / "PACE.top"
    shutil.copy2(topology_top, pace_top_dest)
    print(f"Copied PACE topology: {topology_top} -> {pace_top_dest}")
    
    # Copy processed structure as conf.gro to openmm directory
    conf_gro_dest = openmm_dir / "conf.gro"
    shutil.copy2(processed_pdb, conf_gro_dest)
    print(f"Copied structure as conf.gro: {processed_pdb} -> {conf_gro_dest}")
    
    task_type = config.get('task_type', 'IDP')
    optimized_structures = run_openmm_optimization(conf_gro_dest, pace_top_dest, OUTPUT_DIR, task_type)
    
    print("="*60)
    print("OpenMM optimization workflow completed!")
    print(f"Task type: {task_type}")
    print(f"Generated topology: {topology_top}")
    print(f"Input structure: {input_structure}")
    print(f"OpenMM files: {conf_gro_dest}, {pace_top_dest}")
    if optimized_structures:
        print(f"Optimized structures: {optimized_structures}")
    print("="*60)

if __name__ == "__main__":
    run_openmm_workflow()
