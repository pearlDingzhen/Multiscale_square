#!/usr/bin/env python3
"""
Stage 4: Solvent Preparation
Auto-generated script for building and equilibrating explicit solvent system with GROMACS.
"""

import os
import sys
import yaml
import shutil
from pathlib import Path

# Import required modules
from multiscale2.utils import run_command

# ============================================================================
# USER CONFIGURABLE VARIABLES
# ============================================================================
# These variables can be modified by the user if needed

# Input directory from Stage 3 (leave None to auto-detect)
INPUT_DIR = None

# Output directory for solvent preparation
OUTPUT_DIR = "output_solvent"

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
    """Auto-detect input directory from Stage 3 output."""
    openmm_dir = "output_openmm"
    if os.path.exists(openmm_dir):
        return openmm_dir
    
    print(f"Default OpenMM directory not found: {openmm_dir}")
    return None

class SolventPreparator:
    """Handles explicit solvent system preparation and equilibration."""
    
    def __init__(self, config, working_dir):
        """Initialize solvent preparator with configuration."""
        self.working_dir = Path(working_dir)
        self.working_dir.mkdir(exist_ok=True)
        
        # Extract key parameters from config
        self.protein_name = config['protein']['name']
        self.nmol = config['protein']['nmol']
        self.task_type = config.get('task_type', 'IDP')
        
        # Extract GROMACS parameters from config
        gromacs_config = config.get('gromacs_explicit_solvent_build_and_equilibration', {})
        self.num_cpus = gromacs_config.get('num_cpus', 1)
        self.num_mpi_threads = gromacs_config.get('num_mpi_threads', 1)
        self.num_omp_threads = gromacs_config.get('num_omp_threads', self.num_cpus)
        self.gmx_execute = gromacs_config.get('gmx_excute', None)
        self.ions_concentration = gromacs_config.get('ions_concentration', 0.15)
        self.temperature = gromacs_config.get('temperature', 300.0)
        self.pressure = gromacs_config.get('pressure', 1.0)
        self.use_anti_freeze_water = gromacs_config.get('Use_anti_freeze_water', True)
        self.anti_free_water_ratio = gromacs_config.get('Anti_free_water_ratio', 0.1)
        self.water_gro = gromacs_config.get('water_gro', 'cg216water.gro')
        self.constraints = gromacs_config.get('constraints', None)
        
        # Try to import GromacsWrapper
        try:
            import gromacs
            self.HAS_GROMACS_WRAPPER = True
        except ImportError:
            self.HAS_GROMACS_WRAPPER = False
        
        # Determine GROMACS execution method
        if self.gmx_execute:
            self.use_gromacs_wrapper = False
            print(f"Using manual GROMACS execution: {self.gmx_execute}")
        elif self.HAS_GROMACS_WRAPPER:
            self.use_gromacs_wrapper = True
            print("Using GROMACS wrapper for execution")
        else:
            raise RuntimeError(
                "Neither gmx_excute specified in config nor GromacsWrapper available. "
                "Please either specify gmx_excute in config or install GromacsWrapper."
            )
        
        # Set up subdirectories
        self.vacuum_dir = self.working_dir / "vacuum_optimization"
        self.solvent_dir = self.working_dir / "explicit_solvent"
        
        for dir_path in [self.vacuum_dir, self.solvent_dir]:
            dir_path.mkdir(exist_ok=True)

    def _run_gromacs_command(self, command, cwd=None, input_data=None):
        """
        Run GROMACS command using either wrapper or manual execution.
        
        Args:
            command: Command to run (without 'gmx' prefix for wrapper)
            cwd: Working directory
            input_data: Input data for interactive commands
        """
        if self.use_gromacs_wrapper:
            # Store current directory and change to target directory if specified
            original_cwd = os.getcwd()
            try:
                if cwd:
                    os.chdir(cwd)
                
                # Use gromacs wrapper
                if 'grompp' in command:
                    # Parse grompp command
                    parts = command.split()
                    args = {}
                    for i in range(len(parts)):
                        if parts[i] == '-f':
                            args['f'] = parts[i+1]
                        elif parts[i] == '-c':
                            args['c'] = parts[i+1]
                        elif parts[i] == '-p':
                            args['p'] = parts[i+1]
                        elif parts[i] == '-o':
                            args['o'] = parts[i+1]
                        elif parts[i] == '-maxwarn':
                            args['maxwarn'] = int(parts[i+1])
                    
                    import gromacs
                    gromacs.grompp(**args)
                    
                elif 'mdrun' in command:
                    # Parse mdrun command
                    parts = command.split()
                    args = {}
                    for i in range(len(parts)):
                        if parts[i] == '-deffnm':
                            args['deffnm'] = parts[i+1]
                        elif parts[i] == '-ntmpi':
                            args['ntmpi'] = int(parts[i+1])
                        elif parts[i] == '-ntomp':
                            args['ntomp'] = int(parts[i+1])
                        elif parts[i] == '-nt':
                            args['nt'] = int(parts[i+1])
                        elif parts[i] == '-v':
                            args['v'] = True
                    
                    import gromacs
                    gromacs.mdrun(**args)
                    
                elif 'solvate' in command:
                    # Parse solvate command
                    parts = command.split()
                    args = {}
                    for i in range(len(parts)):
                        if parts[i] == '-cp':
                            args['cp'] = parts[i+1]
                        elif parts[i] == '-cs':
                            args['cs'] = parts[i+1]
                        elif parts[i] == '-o':
                            args['o'] = parts[i+1]
                        elif parts[i] == '-p':
                            args['p'] = parts[i+1]
                    
                    import gromacs
                    gromacs.solvate(**args)
                    
                elif 'genion' in command:
                    # Parse genion command
                    parts = command.split()
                    args = {}
                    for i in range(len(parts)):
                        if parts[i] == '-s':
                            args['s'] = parts[i+1]
                        elif parts[i] == '-o':
                            args['o'] = parts[i+1]
                        elif parts[i] == '-p':
                            args['p'] = parts[i+1]
                        elif parts[i] == '-conc':
                            args['conc'] = float(parts[i+1])
                        elif parts[i] == '-neutral':
                            args['neutral'] = True
                    
                    import gromacs
                    if input_data:
                        gromacs.genion(input=input_data, **args)
                    else:
                        gromacs.genion(**args)
                else:
                    # Fallback to run_command for unsupported commands
                    full_command = f"gmx {command}"
                    run_command(full_command, cwd=cwd)
            finally:
                # Always restore original directory
                os.chdir(original_cwd)
        else:
            # Use manual execution
            full_command = f"{self.gmx_execute} {command}"
            if input_data:
                # Handle interactive commands with input data
                import subprocess
                import shlex
                import logging
                
                logging.info(f"Executing command in '{cwd}':\n  $ {full_command}")
                logging.info(f"With input: {input_data}")
                
                try:
                    result = subprocess.run(
                        shlex.split(full_command),
                        cwd=cwd,
                        check=True,
                        text=True,
                        capture_output=True,
                        input=input_data
                    )
                    if result.stdout:
                        logging.info(f"STDOUT:\n{result.stdout}")
                    if result.stderr:
                        logging.warning(f"STDERR:\n{result.stderr}")
                except subprocess.CalledProcessError as e:
                    logging.error(f"Command failed with exit code {e.returncode}")
                    logging.error(f"STDOUT:\n{e.stdout}")
                    logging.error(f"STDERR:\n{e.stderr}")
                    raise
            else:
                run_command(full_command, cwd=cwd)

    def _write_mdp(self, filename, content, target_dir):
        """Write MDP file content to specified directory."""
        mdp_path = target_dir / filename
        with open(mdp_path, 'w') as f:
            f.write(content)
        print(f"✓ Created MDP file: {mdp_path}")
    
    def _get_constraints_string(self):
        """Get constraints string based on configuration."""
        if self.constraints is None or self.constraints == 'null':
            return 'none'
        elif self.constraints == 'hbonds':
            return 'h-bonds'
        else:
            # Default to none for unknown values
            return 'none'

    def _create_vacuum_em_mdps(self):
        """Create MDP files for vacuum energy minimization."""
        
        constraints_str = self._get_constraints_string()
        
        # Vacuum steep descent minimization
        em_steep_mdp = f"""define              =  
constraints         =  {constraints_str}
integrator          =  steep
nsteps              =  50000
lincs_iter          =  8
emtol               =  500
emstep              =  0.01

cutoff-scheme = Verlet
nstlist             =  10
ns_type             =  grid
rlist               =  1.2
coulombtype         =  reaction-field
rcoulomb = 1.1
epsilon_r           =  15
vdwtype = Cut-off
vdw_modifier=Force-switch
rvdw                =  1.1
rvdw_switch         =  0.9

nstcomm             =  100
table-extension     =  10
Tcoupl              =  no
Pcoupl              =  no
gen_vel             =  no
"""
        
        # Vacuum conjugate gradient minimization  
        em_cg_mdp = f"""define              =  
constraints         =  {constraints_str}
integrator          =  cg
nsteps              =  25000
lincs_iter          =  8
emtol               =  100
emstep              =  0.01

cutoff-scheme = Verlet
nstlist             =  10
ns_type             =  grid
rlist               =  1.2
coulombtype         =  reaction-field
rcoulomb = 1.1
epsilon_r           =  15
vdwtype = Cut-off
vdw_modifier=Force-switch
rvdw                =  1.1
rvdw_switch         =  0.9

nstcomm             =  100
table-extension     =  10
Tcoupl              =  no
Pcoupl              =  no
gen_vel             =  no
"""
        
        self._write_mdp('em_steep.mdp', em_steep_mdp, self.vacuum_dir)
        self._write_mdp('em_cg.mdp', em_cg_mdp, self.vacuum_dir)

    def _create_solvent_em_mdps(self):
        """Create MDP files for explicit solvent energy minimization."""
        
        constraints_str = self._get_constraints_string()
        
        # Solvent steep descent minimization
        em_steep_mdp = f"""define              =  
constraints         =  {constraints_str}
integrator          =  steep
nsteps              =  50000
lincs_iter          =  8
emtol               =  2000
emstep              =  0.01

cutoff-scheme = Verlet
nstlist             =  10
ns_type             =  grid
rlist               =  1.2
coulombtype         =  reaction-field
rcoulomb = 1.1
epsilon_r           =  15
vdwtype = Cut-off
vdw_modifier=Force-switch
rvdw                =  1.1
rvdw_switch         =  0.9

nstcomm             =  100
table-extension     =  10
Tcoupl              =  no
Pcoupl              =  no
gen_vel             =  no
"""
        
        # Solvent conjugate gradient minimization  
        em_cg_mdp = f"""define              =  
constraints         =  {constraints_str}
integrator          =  cg
nsteps              =  25000
lincs_iter          =  8
emtol               =  200
emstep              =  0.01

cutoff-scheme = Verlet
nstlist             =  10
ns_type             =  grid
rlist               =  1.2
coulombtype         =  reaction-field
rcoulomb = 1.1
epsilon_r           =  15
vdwtype = Cut-off
vdw_modifier=Force-switch
rvdw                =  1.1
rvdw_switch         =  0.9

nstcomm             =  100
table-extension     =  10
Tcoupl              =  no
Pcoupl              =  no
gen_vel             =  no
"""
        
        self._write_mdp('em_steep.mdp', em_steep_mdp, self.solvent_dir)
        self._write_mdp('em_cg.mdp', em_cg_mdp, self.solvent_dir)

    def step1_vacuum_optimization(self, openmm_output_dir):
        """
        Step 1: Optimize protein in vacuum.
        
        Args:
            openmm_output_dir: Path to the OpenMM output directory
            
        Returns:
            tuple: (optimized_gro_path, topology_path) - paths to optimized files
        """
        print("="*60)
        print("Step 1: Vacuum Protein Optimization")
        print("="*60)
        
        # Get OpenMM outputs
        openmm_dir = Path(openmm_output_dir).resolve()
        conf_gro_path = openmm_dir / "conf.gro"
        pace_top_path = openmm_dir / "PACE.top"
        
        if not conf_gro_path.exists():
            raise FileNotFoundError(f"conf.gro not found at: {conf_gro_path}")
        if not pace_top_path.exists():
            raise FileNotFoundError(f"PACE.top not found at: {pace_top_path}")
        
        print(f"✓ Found OpenMM outputs:")
        print(f"  - {conf_gro_path}")
        print(f"  - {pace_top_path}")
        
        # Copy files to vacuum optimization directory
        vacuum_conf_gro = self.vacuum_dir / "conf.gro"
        vacuum_pace_top = self.vacuum_dir / "PACE.top"
        
        shutil.copy2(conf_gro_path, vacuum_conf_gro)
        shutil.copy2(pace_top_path, vacuum_pace_top)
        
        print(f"✓ Copied files to vacuum optimization directory")
        
        # Create MDP files for vacuum optimization
        self._create_vacuum_em_mdps()
        
        # Run vacuum energy minimization
        print("Running vacuum energy minimization...")
        print("  - Steep descent minimization (emtol=500)...")
        
        # Steep descent minimization
        self._run_gromacs_command(
            f"grompp -f em_steep.mdp -c conf.gro -p PACE.top -o em_steep.tpr -maxwarn 2",
            cwd=self.vacuum_dir
        )
        
        self._run_gromacs_command(
            f"mdrun -v -deffnm em_steep -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
            cwd=self.vacuum_dir
        )
        
        print("  - Conjugate gradient minimization (emtol=100)...")
        
        # Conjugate gradient minimization
        self._run_gromacs_command(
            f"grompp -f em_cg.mdp -c em_steep.gro -p PACE.top -o em_cg.tpr -maxwarn 2",
            cwd=self.vacuum_dir
        )
        
        self._run_gromacs_command(
            f"mdrun -v -deffnm em_cg -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
            cwd=self.vacuum_dir
        )
        
        # Return paths to optimized files
        optimized_gro = self.vacuum_dir / "em_cg.gro"
        optimized_top = self.vacuum_dir / "PACE.top"
        
        print(f"✓ Vacuum optimization completed!")
        print(f"✓ Optimized structure: {optimized_gro}")
        print(f"✓ Topology file: {optimized_top}")
        
        return optimized_gro, optimized_top

    def step2_build_solvent_system(self, optimized_gro_path, topology_path):
        """
        Step 2: Build explicit solvent system with anti-freeze water if needed.
        
        Args:
            optimized_gro_path: Path to the vacuum-optimized structure
            topology_path: Path to the topology file
            
        Returns:
            tuple: (solvated_gro_path, solvated_top_path) - paths to solvated files
        """
        print("="*60)
        print("Step 2: Build Explicit Solvent System")
        print("="*60)
        
        # Copy optimized files to solvent directory
        solvent_gro = self.solvent_dir / "system.gro"
        solvent_top = self.solvent_dir / "topol.top"
        
        shutil.copy2(optimized_gro_path, solvent_gro)
        shutil.copy2(topology_path, solvent_top)
        
        print(f"✓ Copied optimized files to solvent directory")
        
        # Create MDP files for solvent system
        self._create_solvent_em_mdps()
        
        # Solvate the system
        print(f"Solvating the system with {self.water_gro}...")
        self._run_gromacs_command(
            f"solvate -cp system.gro -cs {self.water_gro} -o solvated.gro -p topol.top",
            cwd=self.solvent_dir
        )
        
        # Add anti-freeze water if enabled
        if self.use_anti_freeze_water:
            print(f"Adding anti-freeze water (ratio: {self.anti_free_water_ratio})...")
            self._add_anti_freeze_water()
        
        # Add ions
        print(f"Adding ions (concentration: {self.ions_concentration} M)...")
        self._run_gromacs_command(
            f"grompp -f em_steep.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2",
            cwd=self.solvent_dir
        )
        
        # Generate ions
        self._run_gromacs_command(
            f"genion -s ions.tpr -o system.gro -p topol.top -conc {self.ions_concentration} -neutral",
            cwd=self.solvent_dir,
            input_data="SOL"
        )
        
        # Final solvated files
        solvated_gro = self.solvent_dir / "system.gro"
        solvated_top = self.solvent_dir / "topol.top"
        
        print(f"✓ Solvent system built successfully!")
        print(f"✓ Solvated structure: {solvated_gro}")
        print(f"✓ Topology file: {solvated_top}")
        
        return solvated_gro, solvated_top

    def _add_anti_freeze_water(self):
        """Add anti-freeze water molecules to the system."""
        # This is a placeholder - actual implementation would depend on
        # the specific anti-freeze water model being used
        print(f"  - Adding {self.anti_free_water_ratio * 100}% anti-freeze water...")
        # Implementation would go here
        pass

    def step3_equilibration(self, solvated_gro_path, solvated_top_path):
        """
        Step 3: Equilibrate the solvated system.
        
        Args:
            solvated_gro_path: Path to the solvated structure
            solvated_top_path: Path to the solvated topology
            
        Returns:
            tuple: (equilibrated_gro_path, equilibrated_top_path) - paths to equilibrated files
        """
        print("="*60)
        print("Step 3: System Equilibration")
        print("="*60)
        
        # Create equilibration MDP files
        self._create_equilibration_mdps()
        
        # Energy minimization
        print("Running energy minimization...")
        self._run_gromacs_command(
            f"grompp -f em_steep.mdp -c system.gro -p topol.top -o em_steep.tpr -maxwarn 2",
            cwd=self.solvent_dir
        )
        
        self._run_gromacs_command(
            f"mdrun -v -deffnm em_steep -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
            cwd=self.solvent_dir
        )
        
        self._run_gromacs_command(
            f"grompp -f em_cg.mdp -c em_steep.gro -p topol.top -o em_cg.tpr -maxwarn 2",
            cwd=self.solvent_dir
        )
        
        self._run_gromacs_command(
            f"mdrun -v -deffnm em_cg -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
            cwd=self.solvent_dir
        )
        
        # NVT equilibration
        print("Running NVT equilibration...")
        self._run_gromacs_command(
            f"grompp -f nvt.mdp -c em_cg.gro -p topol.top -o nvt.tpr -maxwarn 2",
            cwd=self.solvent_dir
        )
        
        self._run_gromacs_command(
            f"mdrun -v -deffnm nvt -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
            cwd=self.solvent_dir
        )
        
        # NPT equilibration
        print("Running NPT equilibration...")
        self._run_gromacs_command(
            f"grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 2",
            cwd=self.solvent_dir
        )
        
        self._run_gromacs_command(
            f"mdrun -v -deffnm npt -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
            cwd=self.solvent_dir
        )
        
        # Final equilibrated files
        equilibrated_gro = self.solvent_dir / "npt.gro"
        equilibrated_top = self.solvent_dir / "topol.top"
        
        print(f"✓ Equilibration completed!")
        print(f"✓ Equilibrated structure: {equilibrated_gro}")
        print(f"✓ Topology file: {equilibrated_top}")
        
        return equilibrated_gro, equilibrated_top

    def _create_equilibration_mdps(self):
        """Create MDP files for equilibration."""
        constraints_str = self._get_constraints_string()
        
        # NVT equilibration
        nvt_mdp = f"""define              = -DPOSRES
constraints         = {constraints_str}
integrator          = md
nsteps              = 50000
dt                  = 0.002

cutoff-scheme = Verlet
nstlist             = 10
ns_type             = grid
rlist               = 1.2
coulombtype         = reaction-field
rcoulomb = 1.1
epsilon_r           = 15
vdwtype = Cut-off
vdw_modifier=Force-switch
rvdw                = 1.1
rvdw_switch         = 0.9

tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1 0.1
ref_t               = {self.temperature} {self.temperature}

pcoupl              = no
gen_vel             = yes
gen_temp            = {self.temperature}
gen_seed            = -1
"""
        
        # NPT equilibration
        npt_mdp = f"""define              = -DPOSRES
constraints         = {constraints_str}
integrator          = md
nsteps              = 50000
dt                  = 0.002

cutoff-scheme = Verlet
nstlist             = 10
ns_type             = grid
rlist               = 1.2
coulombtype         = reaction-field
rcoulomb = 1.1
epsilon_r           = 15
vdwtype = Cut-off
vdw_modifier=Force-switch
rvdw                = 1.1
rvdw_switch         = 0.9

tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1 0.1
ref_t               = {self.temperature} {self.temperature}

pcoupl              = Parrinello-Rahman
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = {self.pressure}
compressibility     = 4.5e-5

gen_vel             = no
"""
        
        self._write_mdp('nvt.mdp', nvt_mdp, self.solvent_dir)
        self._write_mdp('npt.mdp', npt_mdp, self.solvent_dir)

    def run_solvent_workflow(self, openmm_output_dir):
        """
        Run the complete solvent preparation workflow.
        
        Args:
            openmm_output_dir: Path to the OpenMM output directory
            
        Returns:
            dict: Results of the solvent preparation workflow
        """
        print("="*80)
        print("Starting Solvent Preparation Workflow for", self.protein_name)
        print(f"Number of molecules: {self.nmol}")
        print(f"Task type: {self.task_type}")
        print(f"Temperature: {self.temperature} K")
        print(f"CPUs: {self.num_cpus}")
        print(f"MPI threads: {self.num_mpi_threads}")
        print(f"OpenMP threads: {self.num_omp_threads}")
        print(f"Anti-freeze water: {self.use_anti_freeze_water} (ratio: {self.anti_free_water_ratio})")
        print("="*80)
        
        try:
            # Step 1: Vacuum optimization
            optimized_gro, optimized_top = self.step1_vacuum_optimization(openmm_output_dir)
            
            # Step 2: Build solvent system
            solvated_gro, solvated_top = self.step2_build_solvent_system(optimized_gro, optimized_top)
            
            # Step 3: Equilibration
            equilibrated_gro, equilibrated_top = self.step3_equilibration(solvated_gro, solvated_top)
            
            print("="*80)
            print("Solvent preparation workflow completed successfully!")
            print("="*80)
            
            return {
                'status': 'success',
                'vacuum_optimized': {
                    'gro': optimized_gro,
                    'top': optimized_top
                },
                'solvated': {
                    'gro': solvated_gro,
                    'top': solvated_top
                },
                'equilibrated': {
                    'gro': equilibrated_gro,
                    'top': equilibrated_top
                }
            }
            
        except Exception as e:
            print(f"✗ Error in solvent preparation: {e}")
            import traceback
            traceback.print_exc()
            return {
                'status': 'error',
                'error': str(e)
            }

def run_solvent_workflow():
    """Execute the complete solvent preparation workflow."""
    print("="*60)
    print("Stage 4: Solvent Preparation Workflow")
    print("="*60)
    
    # Load configuration
    config = load_config()
    
    # Determine input directory
    input_dir = INPUT_DIR or auto_detect_input_dir()
    if not input_dir:
        print("Error: Could not determine input directory from Stage 3")
        print("Please set INPUT_DIR variable or ensure Stage 3 completed successfully")
        sys.exit(1)
    
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Initialize solvent preparator
    preparator = SolventPreparator(config, OUTPUT_DIR)
    
    # Run the workflow
    results = preparator.run_solvent_workflow(input_dir)
    
    if results['status'] == 'success':
        print("\n✓ Solvent preparation completed successfully!")
        print(f"Final equilibrated structure: {results['equilibrated']['gro']}")
        print(f"Final topology: {results['equilibrated']['top']}")
    else:
        print(f"\n✗ Solvent preparation failed: {results.get('error', 'Unknown error')}")
        sys.exit(1)

if __name__ == "__main__":
    run_solvent_workflow()
