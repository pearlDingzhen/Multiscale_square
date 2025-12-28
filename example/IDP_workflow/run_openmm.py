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
        # Step 1: Minimize without nonbonded interactions
        ################################################################################
        
        print("Step 1: Minimizing with bond-only interactions...")
        
        system = top.create_system(nonbonded_cutoff=1.1 * unit.nanometer)
        system.removeForce(0)  # Remove nonbonded force
        
        integrator = LangevinIntegrator(ref_t * unit.kelvin,
                                      10.0 / unit.picosecond,
                                      time_step * unit.picosecond)
        integrator.setRandomNumberSeed(0)
        
        # Add CA restraint only for IDP tasks
        if task_type == 'IDP':
            print("Adding CA restraints for IDP task...")
            restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
            system.addForce(restraint)
            restraint.addGlobalParameter('k', 100.0*unit.kilojoules_per_mole/unit.nanometer)
            restraint.addPerParticleParameter('x0')
            restraint.addPerParticleParameter('y0')
            restraint.addPerParticleParameter('z0')
            
            for atom in top.topology.atoms():
                if atom.name == 'CA':
                    restraint.addParticle(atom.index, gro_position[atom.index])
        
        simulation = Simulation(top.topology, system, integrator, platform)
        simulation.context.setPositions(conf.getPositions())
        
        # Energy before minimization
        energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        print(f'Energy (without nonbond) before minimization: {energies}')
        
        # Minimize
        simulation.minimizeEnergy(maxIterations=5000, tolerance=100)
        
        # Energy after minimization
        energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        state_bond = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        PDBFile.writeFile(simulation.topology, state_bond.getPositions(), open('optibond.pdb', 'w'))
        
        print(f'Energy (without nonbond) after minimization: {energies}')
        print('Bond-only minimized structure saved as optibond.pdb')
        
        ################################################################################
        # Step 2: Minimize with full interactions
        ################################################################################
        
        print("\nStep 2: Minimizing with full interactions...")
        
        system1 = top.create_system(nonbonded_cutoff=1.1 * unit.nanometer)
        
        # Add CA restraint only for IDP tasks
        if task_type == 'IDP':
            restraint1 = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
            system1.addForce(restraint1)
            restraint1.addGlobalParameter('k', 100*unit.kilojoules_per_mole/unit.nanometer)
            restraint1.addPerParticleParameter('x0')
            restraint1.addPerParticleParameter('y0')
            restraint1.addPerParticleParameter('z0')
            
            for atom in top.topology.atoms():
                if atom.name == 'CA':
                    restraint1.addParticle(atom.index, gro_position[atom.index])
        
        integrator1 = LangevinIntegrator(ref_t * unit.kelvin,
                                       10.0 / unit.picosecond,
                                       time_step * unit.picosecond)
        integrator1.setRandomNumberSeed(0)
        
        simulation1 = Simulation(top.topology, system1, integrator1, platform, properties)
        simulation1.context.setPositions(state_bond.getPositions())
        
        # Energy before minimization
        state = simulation1.context.getState(getEnergy=True, getForces=False)
        print(f'Energy (with nonbond) before minimization: {state.getPotentialEnergy()}')
        
        # Minimize
        simulation1.minimizeEnergy(maxIterations=5000, tolerance=200)
        
        # Final state
        state_nonbond = simulation1.context.getState(getPositions=True, enforcePeriodicBox=True)
        PDBFile.writeFile(simulation1.topology, state_nonbond.getPositions(), open('optinonbond.pdb', 'w'))
        
        state = simulation1.context.getState(getEnergy=True, getForces=False)
        print(f'Energy (with nonbond) after minimization: {state.getPotentialEnergy()}')
        
        print('\033[32mMinimization successful! Final structure saved as optinonbond.pdb\033[0m')
        
        # Check for output files
        output_files = ["optibond.pdb", "optinonbond.pdb"]
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
