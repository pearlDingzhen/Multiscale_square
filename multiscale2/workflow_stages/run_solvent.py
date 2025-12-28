#!/usr/bin/env python3
"""
Stage 4: Solvent Preparation
Auto-generated script for building and equilibrating explicit solvent system with GROMACS.

REQUIRED INPUT FILES FROM STAGE 3 (OpenMM optimization):
The script expects to find the following files in the OpenMM output directory:
- opti2.pdb: Optimized structure file from OpenMM optimization (final structure)
- PACE.top: Topology file from PACE topology generation

SEARCH LOCATIONS (in order):
1. output_openmm/openmm/opti2.pdb and output_openmm/openmm/PACE.top
2. output_openmm/opti2.pdb and output_openmm/PACE.top  
3. output_openmm/topology/opti2.pdb and output_openmm/topology/PACE.top
4. output_openmm/system.pdb (as fallback for structure)

If files are not found, ensure Stage 3 (OpenMM optimization) completed successfully.
"""

import os
import sys
import yaml
import shutil
import re
import subprocess
import time
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
    """Auto-detect input directory from Stage 3 output by finding opti2.pdb."""
    # Try to find opti2.pdb in common locations
    possible_dirs = [
        "output_openmm/openmm",  # Most common location
        "output_openmm",          # Alternative location
        ".",                      # Current directory
    ]
    
    for openmm_dir in possible_dirs:
        opti2_pdb = os.path.join(openmm_dir, "opti2.pdb")
        if os.path.exists(opti2_pdb):
            print(f"✓ Found opti2.pdb in: {openmm_dir}")
            return openmm_dir
    
    print(f"❌ opti2.pdb not found!")
    print(f"   Searched in: {possible_dirs}")
    print(f"   Please ensure Stage 3 (OpenMM optimization) completed successfully")
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
        
        # Equilibration parameters
        self.nvt_steps = gromacs_config.get('NVT_steps', 1000000)
        self.npt_steps = gromacs_config.get('NPT_steps', 1000000)
        self.dt = gromacs_config.get('dt', 0.002)
        self.use_gpu = gromacs_config.get('use_gpu', False)
        
        # Get topology type from CG CALVADOS config
        cg_config = config.get('cg_calvados', {})
        self.topol = cg_config.get('topol', 'grid')
        
        # EM define configuration (for RINGREFINE and other preprocessor definitions)
        # Default: use RINGREFINE_HEAVY for first stage, none for second stage
        self.em_define_primary = '-DRINGREFINE_HEAVY -DRINGREFINE'
        self.em_define_secondary = ''  # Empty = no define for second stage
        
        # EM retry configuration
        self.em_max_retries = gromacs_config.get('em_max_retries', 5)
        
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

    def _check_em_success(self, log_file):
        """
        Check if energy minimization was successful by parsing the log file.
        
        Args:
            log_file: Path to the GROMACS log file
            
        Returns:
            tuple: (success: bool, error_message: str)
        """
        if not os.path.exists(log_file):
            return False, f"Log file not found: {log_file}"
        
        try:
            with open(log_file, 'r') as f:
                content = f.read()
            
            # Check for convergence indicators
            convergence_patterns = [
                r'converged to Fmax < \d+',
                r'converged to machine precision',
                r'Energy minimization has converged',
                r'Steepest Descents converged',
                r'Polak-Ribiere Conjugate Gradients converged'
            ]
            
            # Check for failure indicators
            failure_patterns = [
                r'but did not reach the requested Fmax',
                r'Energy minimization has stopped',
                r'Maximum number of steps reached',
                r'Fatal error',
                r'Error'
            ]
            
            # Look for convergence
            for pattern in convergence_patterns:
                if re.search(pattern, content, re.IGNORECASE):
                    # Check if it's a true convergence or just partial
                    if 'but did not reach the requested Fmax' in content:
                        return False, "Converged to machine precision but did not reach requested Fmax"
                    return True, "Energy minimization converged successfully"
            
            # Look for explicit failures
            for pattern in failure_patterns:
                if re.search(pattern, content, re.IGNORECASE):
                    return False, f"EM failed: {pattern}"
            
            # If no clear success or failure, check the end of the file
            lines = content.strip().split('\n')
            if lines:
                last_line = lines[-1]
                if 'Finished mdrun' in last_line:
                    return True, "EM completed (no explicit convergence message)"
            
            return False, "No clear success or failure indication found"
            
        except Exception as e:
            return False, f"Error reading log file: {e}"

    def _run_em_with_retry(self, mdp_file, input_gro, topology, output_prefix, cwd=None):
        """
        Run energy minimization with automatic retry on failure.
        
        Args:
            mdp_file: MDP file path
            input_gro: Input structure file
            topology: Topology file
            output_prefix: Output file prefix
            cwd: Working directory
            
        Returns:
            tuple: (success: bool, final_output: str, attempts: int)
        """
        if cwd is None:
            cwd = os.getcwd()
        
        print(f"Starting energy minimization with retry mechanism...")
        print(f"  Input: {input_gro}")
        print(f"  MDP: {mdp_file}")
        print(f"  Topology: {topology}")
        print(f"  Output prefix: {output_prefix}")
        print(f"  Max retries: {self.em_max_retries}")
        print(f"  GROMACS: {self.gmx_execute if self.gmx_execute else 'gromacsWrapper'}")
        print()
        
        for attempt in range(self.em_max_retries + 1):
            if attempt > 0:
                print(f"🔄 Retry attempt {attempt}/{self.em_max_retries}")
                time.sleep(2)  # Brief pause between retries
            
            print(f"📋 Attempt {attempt + 1}: Running energy minimization...")
            
            try:
                # Generate commands
                grompp_cmd, mdrun_cmd = self._get_em_commands(
                    mdp_file, input_gro, topology, output_prefix, cwd
                )
                
                # Run grompp
                print(f"  Running: {' '.join(grompp_cmd)}")
                result = subprocess.run(grompp_cmd, cwd=cwd, capture_output=True, text=True)
                if result.returncode != 0:
                    print(f"  ❌ grompp failed: {result.stderr}")
                    continue
                
                # Run mdrun
                print(f"  Running: {' '.join(mdrun_cmd)}")
                result = subprocess.run(mdrun_cmd, cwd=cwd, capture_output=True, text=True)
                if result.returncode != 0:
                    print(f"  ❌ mdrun failed: {result.stderr}")
                    continue
                
                # Check log file
                log_file = os.path.join(cwd, f"{output_prefix}.log")
                success, message = self._check_em_success(log_file)
                
                if success:
                    print(f"  ✅ Success: {message}")
                    return True, output_prefix, attempt + 1
                else:
                    print(f"  ❌ Failed: {message}")
                    if attempt < self.em_max_retries:
                        print(f"  🔄 Will retry...")
                    else:
                        print(f"  🛑 Maximum retries reached")
                
            except Exception as e:
                print(f"  ❌ Exception during EM: {e}")
                if attempt < self.em_max_retries:
                    print(f"  🔄 Will retry...")
                else:
                    print(f"  🛑 Maximum retries reached")
        
        # All attempts failed
        print(f"\n❌ Energy minimization failed after {self.em_max_retries + 1} attempts")
        print(f"💡 Suggestions:")
        print(f"   1. Check your input structure for clashes or unrealistic geometry")
        print(f"   2. Try using double precision GROMACS (set gmx_excute to 'gmx_d')")
        print(f"   3. Increase emtol in your MDP file")
        print(f"   4. Check your force field parameters")
        
        return False, output_prefix, self.em_max_retries + 1

    def _get_em_commands(self, mdp_file, input_gro, topology, output_prefix, cwd=None):
        """
        Generate the energy minimization commands.
        
        Args:
            mdp_file: MDP file path
            input_gro: Input structure file
            topology: Topology file
            output_prefix: Output file prefix
            cwd: Working directory
            
        Returns:
            tuple: (grompp_cmd, mdrun_cmd)
        """
        # Determine GROMACS executable based on configuration
        if self.gmx_execute:
            gmx_cmd = self.gmx_execute
        else:
            gmx_cmd = "gmx"  # Default for gromacsWrapper
        
        # Use relative paths if cwd is specified
        if cwd:
            mdp_file = os.path.basename(mdp_file)
            input_gro = os.path.basename(input_gro)
            topology = os.path.basename(topology)
        
        # Generate grompp command
        grompp_cmd = [
            gmx_cmd, "grompp",
            "-f", mdp_file,
            "-c", input_gro,
            "-p", topology,
            "-o", f"{output_prefix}.tpr",
            "-maxwarn", "2"
        ]
        
        # Generate mdrun command
        mdrun_cmd = [
            gmx_cmd, "mdrun",
            "-v",
            "-deffnm", output_prefix,
            "-ntmpi", str(self.num_mpi_threads),
            "-ntomp", str(self.num_omp_threads)
        ]
        
        return grompp_cmd, mdrun_cmd

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

    def _get_define_string(self, define_value=None):
        """Get the define string for MDP files.
        
        Args:
            define_value: Optional define value to use. If None, uses em_define_primary.
            
        Returns:
            str: Define line for MDP file, or empty string if no define.
        """
        define = define_value if define_value is not None else self.em_define_primary
        if define:
            return f"define              =  {define}"
        return ""

    def _get_secondary_define_string(self):
        """Get the define string for secondary (no ringrefine) minimization.
        
        Returns:
            str: Empty string (no define for second stage).
        """
        return ""

    def _create_cg_em_fallback_mdp(self, target_dir, constraints_str, define_str=None):
        """Create fallback CG MDP with higher emtol (1000) for retry.
        
        This is used when the primary emtol=500 fails to converge.
        
        Args:
            target_dir: Target directory to write MDP file
            constraints_str: Constraints string
            define_str: Optional define string (e.g., '-DRINGREFINE_HEAVY')
        """
        define_line = self._get_define_string(define_str)
        
        em_cg_fallback_mdp = f"""{define_line}
constraints         =  {constraints_str}
integrator          =  cg
nsteps              =  25000
lincs_iter          =  8
emtol               =  1000
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
        self._write_mdp('em_cg_fallback.mdp', em_cg_fallback_mdp, target_dir)
        return str(target_dir / 'em_cg_fallback.mdp')

    def _run_cg_em_with_fallback(self, primary_mdp, fallback_mdp, input_gro, topology, output_prefix, cwd=None):
        """
        Run conjugate gradient energy minimization with fallback to higher emtol.
        
        Strategy:
        1. First try with primary MDP (emtol=500)
        2. If it fails, retry with fallback MDP (emtol=1000)
        3. Only fail if both attempts fail
        
        Args:
            primary_mdp: Primary MDP file path (emtol=500)
            fallback_mdp: Fallback MDP file path (emtol=1000)
            input_gro: Input structure file
            topology: Topology file
            output_prefix: Output file prefix
            cwd: Working directory
            
        Returns:
            tuple: (success: bool, final_output: str, emtol_used: int)
        """
        if cwd is None:
            cwd = os.getcwd()
        
        print(f"Starting CG energy minimization with fallback mechanism...")
        print(f"  Input: {input_gro}")
        print(f"  Output prefix: {output_prefix}")
        print(f"  GROMACS: {self.gmx_execute if self.gmx_execute else 'gromacsWrapper'}")
        print()
        
        # First attempt: emtol=500
        print("📋 Attempt 1: Running CG with emtol=500...")
        grompp_cmd, mdrun_cmd = self._get_em_commands(
            primary_mdp, input_gro, topology, output_prefix, cwd
        )
        
        # Run grompp
        print(f"  Running: {' '.join(grompp_cmd)}")
        result = subprocess.run(grompp_cmd, cwd=cwd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ⚠️  grompp failed: {result.stderr}")
            print(f"  🔄 Proceeding to fallback...")
        else:
            # Run mdrun
            print(f"  Running: {' '.join(mdrun_cmd)}")
            result = subprocess.run(mdrun_cmd, cwd=cwd, capture_output=True, text=True)
            
            # Check log file
            log_file = os.path.join(cwd, f"{output_prefix}.log")
            success, message = self._check_em_success(log_file)
            
            if success:
                print(f"  ✅ Success (emtol=500): {message}")
                return True, output_prefix, 500
        
        # Second attempt: emtol=1000
        print()
        print("📋 Attempt 2: Running CG with emtol=1000 (fallback)...")
        grompp_cmd, mdrun_cmd = self._get_em_commands(
            fallback_mdp, input_gro, topology, output_prefix, cwd
        )
        
        # Run grompp
        print(f"  Running: {' '.join(grompp_cmd)}")
        result = subprocess.run(grompp_cmd, cwd=cwd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ❌ grompp failed: {result.stderr}")
            print(f"\n❌ CG energy minimization failed after all attempts")
            print(f"💡 Suggestions:")
            print(f"   1. Check your input structure for clashes or unrealistic geometry")
            print(f"   2. Try using double precision GROMACS (set gmx_excute to 'gmx_d')")
            print(f"   3. Check your force field parameters")
            return False, output_prefix, -1
        
        # Run mdrun
        print(f"  Running: {' '.join(mdrun_cmd)}")
        result = subprocess.run(mdrun_cmd, cwd=cwd, capture_output=True, text=True)
        
        # Check log file
        log_file = os.path.join(cwd, f"{output_prefix}.log")
        success, message = self._check_em_success(log_file)
        
        if success:
            print(f"  ✅ Success (emtol=1000): {message}")
            return True, output_prefix, 1000
        else:
            print(f"  ❌ Failed (emtol=1000): {message}")
            print(f"\n❌ CG energy minimization failed after all attempts")
            print(f"💡 Suggestions:")
            print(f"   1. Check your input structure for clashes or unrealistic geometry")
            print(f"   2. Try using double precision GROMACS (set gmx_excute to 'gmx_d')")
            print(f"   3. Check your force field parameters")
            return False, output_prefix, -1

    def _create_vacuum_em_mdps(self, define_str=None):
        """Create MDP files for vacuum energy minimization.
        
        Note: Using RINGREFINE_HEAVY for additional diagonal distance restraints
        to maintain aromatic ring planarity in vacuum.
        
        Args:
            define_str: Optional define string (e.g., '-DRINGREFINE_HEAVY'). 
                       If None, uses em_define from config.
        """
        
        constraints_str = self._get_constraints_string()
        define_line = self._get_define_string(define_str)
        
        # Vacuum steep descent minimization
        em_steep_mdp = f"""{define_line}
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
        
        # Vacuum conjugate gradient minimization
        em_cg_mdp = f"""{define_line}
constraints         =  {constraints_str}
integrator          =  cg
nsteps              =  25000
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
        
        self._write_mdp('em_steep.mdp', em_steep_mdp, self.vacuum_dir)
        self._write_mdp('em_cg.mdp', em_cg_mdp, self.vacuum_dir)

    def _create_solvent_em_mdps(self):
        """Create MDP files for explicit solvent energy minimization.
        
        Note: Using RINGREFINE_HEAVY for additional diagonal distance restraints
        to maintain aromatic ring planarity in solvent environment.
        
        Creates two sets of MDP files:
        1. Primary set (with ringrefine_heavy): em_steep.mdp, em_cg.mdp
        2. Secondary set (without ringrefine_heavy): em_steep2.mdp, em_cg2.mdp
        """
        
        constraints_str = self._get_constraints_string()
        define_line = self._get_define_string()  # With ringrefine_heavy
        define_line_secondary = self._get_secondary_define_string()  # Without ringrefine_heavy
        
        # Primary: Solvent steep descent minimization (with ringrefine_heavy)
        em_steep_mdp = f"""{define_line}
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
        
        # Primary: Solvent conjugate gradient minimization (with ringrefine_heavy)
        em_cg_mdp = f"""{define_line}
constraints         =  {constraints_str}
integrator          =  cg
nsteps              =  25000
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
        
        # Secondary: Solvent steep descent minimization (without ringrefine_heavy)
        em_steep2_mdp = f"""{define_line_secondary}
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
        
        # Secondary: Solvent conjugate gradient minimization (without ringrefine_heavy)
        em_cg2_mdp = f"""{define_line_secondary}
constraints         =  {constraints_str}
integrator          =  cg
nsteps              =  25000
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
        
        # Write primary MDPs
        self._write_mdp('em_steep.mdp', em_steep_mdp, self.solvent_dir)
        self._write_mdp('em_cg.mdp', em_cg_mdp, self.solvent_dir)
        
        # Write secondary MDPs (without ringrefine_heavy)
        self._write_mdp('em_steep2.mdp', em_steep2_mdp, self.solvent_dir)
        self._write_mdp('em_cg2.mdp', em_cg2_mdp, self.solvent_dir)

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
        
        # Get OpenMM outputs with intelligent path finding
        openmm_dir = Path(openmm_output_dir).resolve()
        
        # Try multiple possible locations for required files
        possible_conf_locations = [
            openmm_dir / "openmm" / "opti2.pdb",
            openmm_dir / "opti2.pdb",
            openmm_dir / "topology" / "opti2.pdb",
            openmm_dir / "system.pdb",  # Alternative: use system.pdb from topology generation
        ]
        
        possible_top_locations = [
            openmm_dir / "PACE.top",
            openmm_dir / "openmm" / "PACE.top",
            openmm_dir / "topology" / "PACE.top",
        ]
        
        # Find opti2.pdb or alternative structure file
        opti2_pdb_path = None
        for path in possible_conf_locations:
            if path.exists():
                opti2_pdb_path = path
                break
        
        # Find PACE.top file
        pace_top_path = None
        for path in possible_top_locations:
            if path.exists():
                pace_top_path = path
                break
        
        # Check if we found the required files
        if not opti2_pdb_path:
            print(f"❌ Required structure file not found!")
            print(f"   Searched locations: {[str(p) for p in possible_conf_locations]}")
            print(f"   Required files: opti2.pdb (or system.pdb as fallback)")
            raise FileNotFoundError(f"Structure file not found in any expected location")
        
        if not pace_top_path:
            print(f"❌ Required topology file not found!")
            print(f"   Searched locations: {[str(p) for p in possible_top_locations]}")
            print(f"   Required files: PACE.top")
            raise FileNotFoundError(f"PACE.top not found in any expected location")
        
        print(f"✓ Found OpenMM outputs:")
        print(f"  - {opti2_pdb_path}")
        print(f"  - {pace_top_path}")
        
        # Copy files to vacuum optimization directory
        vacuum_opti2_pdb = self.vacuum_dir / "opti2.pdb"
        vacuum_pace_top = self.vacuum_dir / "PACE.top"
        
        shutil.copy2(opti2_pdb_path, vacuum_opti2_pdb)
        shutil.copy2(pace_top_path, vacuum_pace_top)
        
        print(f"✓ Copied files to vacuum optimization directory")
        
        # Create MDP files for vacuum optimization (with ringrefine_heavy)
        self._create_vacuum_em_mdps(define_str=self.em_define_primary)
        
        # Run vacuum energy minimization with retry mechanism
        print("Running vacuum energy minimization...")
        print("  - Steep descent minimization (emtol=2000)...")
        
        # Steep descent minimization with retry
        success, output, attempts = self._run_em_with_retry(
            mdp_file=str(self.vacuum_dir / "em_steep.mdp"),
            input_gro=str(self.vacuum_dir / "opti2.pdb"),
            topology=str(self.vacuum_dir / "PACE.top"),
            output_prefix="em_steep",
            cwd=str(self.vacuum_dir)
        )
        
        if not success:
            raise RuntimeError(f"Vacuum steep descent minimization failed after {attempts} attempts")
        
        print("  - Conjugate gradient minimization (emtol=500 → fallback 1000)...")
        
        # Get constraints string for fallback MDP
        constraints_str = self._get_constraints_string()
        
        # Create fallback MDP for CG (with ringrefine_heavy)
        fallback_mdp_path = self._create_cg_em_fallback_mdp(
            self.vacuum_dir, constraints_str, define_str=self.em_define_primary
        )
        
        # Run CG with fallback mechanism
        success, output, emtol_used = self._run_cg_em_with_fallback(
            primary_mdp=str(self.vacuum_dir / "em_cg.mdp"),
            fallback_mdp=fallback_mdp_path,
            input_gro=str(self.vacuum_dir / "em_steep.gro"),
            topology=str(self.vacuum_dir / "PACE.top"),
            output_prefix="em_cg",
            cwd=str(self.vacuum_dir)
        )
        
        if not success:
            raise RuntimeError("Vacuum conjugate gradient minimization failed (tried emtol=500 and 1000)")
        
        print(f"  ✓ CG converged with emtol={emtol_used}")
        
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
        
        # Energy minimization - Stage 1: with RINGREFINE_HEAVY
        print("Running energy minimization (Stage 1: with RINGREFINE_HEAVY)...")
        print("  - Steep descent minimization (emtol=2000)...")
        
        # Steep descent minimization with retry
        success, output, attempts = self._run_em_with_retry(
            mdp_file=str(self.solvent_dir / "em_steep.mdp"),
            input_gro=str(self.solvent_dir / "system.gro"),
            topology=str(self.solvent_dir / "topol.top"),
            output_prefix="em_steep",
            cwd=str(self.solvent_dir)
        )
        
        if not success:
            raise RuntimeError(f"Solvent steep descent minimization (Stage 1) failed after {attempts} attempts")
        
        print("  - Conjugate gradient minimization (emtol=500 → fallback 1000)...")
        
        # Get constraints string for fallback MDP
        constraints_str = self._get_constraints_string()
        
        # Create fallback MDP for CG (with ringrefine_heavy)
        fallback_mdp_path = self._create_cg_em_fallback_mdp(
            self.solvent_dir, constraints_str, define_str=self.em_define_primary
        )
        
        # Run CG with fallback mechanism
        success, output, emtol_used = self._run_cg_em_with_fallback(
            primary_mdp=str(self.solvent_dir / "em_cg.mdp"),
            fallback_mdp=fallback_mdp_path,
            input_gro=str(self.solvent_dir / "em_steep.gro"),
            topology=str(self.solvent_dir / "topol.top"),
            output_prefix="em_cg",
            cwd=str(self.solvent_dir)
        )
        
        if not success:
            raise RuntimeError("Solvent conjugate gradient minimization (Stage 1) failed (tried emtol=500 and 1000)")
        
        print(f"  ✓ Stage 1 CG converged with emtol={emtol_used}")
        
        # Energy minimization - Stage 2: without RINGREFINE_HEAVY
        print("Running energy minimization (Stage 2: without RINGREFINE_HEAVY)...")
        print("  - Steep descent minimization (emtol=2000)...")
        
        # Steep descent minimization with retry (without ringrefine_heavy)
        success, output, attempts = self._run_em_with_retry(
            mdp_file=str(self.solvent_dir / "em_steep2.mdp"),
            input_gro=str(self.solvent_dir / "em_cg.gro"),
            topology=str(self.solvent_dir / "topol.top"),
            output_prefix="em_steep2",
            cwd=str(self.solvent_dir)
        )
        
        if not success:
            raise RuntimeError(f"Solvent steep descent minimization (Stage 2) failed after {attempts} attempts")
        
        print("  - Conjugate gradient minimization (emtol=500 → fallback 1000)...")
        
        # Create fallback MDP for CG (without ringrefine_heavy)
        fallback_mdp_path2 = self._create_cg_em_fallback_mdp(
            self.solvent_dir, constraints_str, define_str=self.em_define_secondary
        )
        
        # Run CG with fallback mechanism (without ringrefine_heavy)
        success, output, emtol_used = self._run_cg_em_with_fallback(
            primary_mdp=str(self.solvent_dir / "em_cg2.mdp"),
            fallback_mdp=fallback_mdp_path2,
            input_gro=str(self.solvent_dir / "em_steep2.gro"),
            topology=str(self.solvent_dir / "topol.top"),
            output_prefix="em_cg2",
            cwd=str(self.solvent_dir)
        )
        
        if not success:
            raise RuntimeError("Solvent conjugate gradient minimization (Stage 2) failed (tried emtol=500 and 1000)")
        
        print(f"  ✓ Stage 2 CG converged with emtol={emtol_used}")
        
        # Final minimized files (use Stage 2 results)
        solvated_gro = self.solvent_dir / "em_cg2.gro"
        solvated_top = self.solvent_dir / "topol.top"
        
        print(f"✓ Solvent system built and minimized successfully!")
        print(f"  Stage 1: with RINGREFINE_HEAVY → em_cg.gro")
        print(f"  Stage 2: without RINGREFINE_HEAVY → em_cg2.gro")
        print(f"✓ Final minimized structure: {solvated_gro}")
        print(f"✓ Topology file: {solvated_top}")
        
        return solvated_gro, solvated_top

    def _add_anti_freeze_water(self):
        """Add anti-freeze water molecules to the system."""
        print(f"  - Adding {self.anti_free_water_ratio * 100}% anti-freeze water...")
        
        # Read the current topology file
        topol_path = self.solvent_dir / "topol.top"
        alllines = []
        
        with open(topol_path, 'r') as f:
            molecule_flag = False
            for line in f:
                if not molecule_flag:
                    alllines.append(line)
                    if line.startswith('[ molecules'):
                        molecule_flag = True
                else:
                    if line.startswith(';') or len(line.split()) < 2:
                        alllines.append(line)
                    else:
                        srl = line.split()
                        if not srl[0] == 'SOL':
                            alllines.append(line)
                        else:
                            SOL_number = int(srl[1])
                            ASOL_number = int(self.anti_free_water_ratio * SOL_number)
                            alllines.append('ASOL {}\n'.format(ASOL_number))
                            alllines.append('SOL  {}\n'.format(SOL_number - ASOL_number))
        
        # Write the modified topology file
        with open(topol_path, 'w') as f:
            for line in alllines:
                f.write(line)
        
        print(f"    - Converted {ASOL_number} SOL molecules to ASOL (anti-freeze water)")
        print(f"    - Remaining SOL molecules: {SOL_number - ASOL_number}")

    def step3_equilibration(self, solvated_gro_path, solvated_top_path):
        """
        Step 3: Equilibrate the solvated system with warmup stages.
        
        Args:
            solvated_gro_path: Path to the solvated structure
            solvated_top_path: Path to the solvated topology
            
        Returns:
            tuple: (equilibrated_gro_path, equilibrated_top_path) - paths to equilibrated files
        """
        print("="*60)
        print("Step 3: System Equilibration with Warmup Stages")
        print("="*60)
        
        # Create equilibration MDP files and manual scripts
        self._create_equilibration_mdps()
        self._create_manual_run_scripts()
        
        # Default behavior: Skip automatic equilibration, user runs manually
        print("⚠️  Skipping automatic equilibration (NVT/NPT)")
        print("⚠️  Energy minimization was completed in Step 2")
        print("⚠️  User must manually run equilibration steps using the generated scripts")
        print("⚠️  Manual run scripts have been created:")
        print(f"    - {self.solvent_dir}/run_nvt.sh")
        print(f"    - {self.solvent_dir}/run_npt.sh")
        print("⚠️  Run these scripts in order: first run_nvt.sh, then run_npt.sh")
        
        # Return the minimized structure from step 2
        print(f"✓ Returning minimized structure from Step 2")
        print("⚠️  Manual equilibration required before production run")
        
        return solvated_gro_path, solvated_top_path

    def _get_equilibration_mdrun_command(self, output_prefix):
        """
        Generate mdrun command for equilibration with optional GPU support.
        
        Args:
            output_prefix: Output file prefix
            
        Returns:
            str: Complete mdrun command string
        """
        gmx_cmd = self.gmx_execute if self.gmx_execute else "gmx"
        
        # Base mdrun command
        cmd_parts = [
            f"{gmx_cmd} mdrun",
            "-v",
            f"-deffnm {output_prefix}",
            f"-ntmpi {self.num_mpi_threads}",
            f"-ntomp {self.num_omp_threads}"
        ]
        
        # Add GPU support if enabled
        if self.use_gpu:
            cmd_parts.extend(["-bonded gpu", "-update gpu"])
        
        return " ".join(cmd_parts)

    def _create_manual_run_scripts(self):
        """Create manual run scripts for NVT and NPT equilibration."""
        print("Creating manual run scripts...")
        
        # Create NVT run script
        nvt_script = self._create_nvt_script()
        nvt_script_path = self.solvent_dir / "run_nvt.sh"
        with open(nvt_script_path, 'w') as f:
            f.write(nvt_script)
        os.chmod(nvt_script_path, 0o755)  # Make executable
        print(f"✓ Created NVT run script: {nvt_script_path}")
        
        # Create NPT run script
        npt_script = self._create_npt_script()
        npt_script_path = self.solvent_dir / "run_npt.sh"
        with open(npt_script_path, 'w') as f:
            f.write(npt_script)
        os.chmod(npt_script_path, 0o755)  # Make executable
        print(f"✓ Created NPT run script: {npt_script_path}")

    def _create_nvt_script(self):
        """Create NVT manual run script."""
        gmx_cmd = self.gmx_execute if self.gmx_execute else "gmx"
        
        # Build mdrun command with optional GPU support
        mdrun_base = f"{gmx_cmd} mdrun -v -deffnm {{output_prefix}} -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}"
        if self.use_gpu:
            mdrun_base += " -bonded gpu -update gpu"
        
        script = f"""#!/bin/bash
# Manual NVT equilibration script
# Generated by run_solvent.py
# 
# This script runs NVT equilibration with warmup stages:
# 1. NVT Warmup1: dt={self.dt/4.0} ps, steps={int(self.nvt_steps * 0.04)}
# 2. NVT Warmup2: dt={self.dt/2.0} ps, steps={int(self.nvt_steps * 0.02)}
# 3. NVT Production: dt={self.dt} ps, steps={self.nvt_steps}

set -e  # Exit on any error

echo "=========================================="
echo "Starting NVT Equilibration with Warmup"
echo "=========================================="
echo "GROMACS command: {gmx_cmd}"
echo "MPI threads: {self.num_mpi_threads}"
echo "OpenMP threads: {self.num_omp_threads}"
echo "Temperature: {self.temperature} K"
echo "GPU acceleration: {'Enabled' if self.use_gpu else 'Disabled'}"
echo ""

# Check if input file exists
if [ ! -f "em_cg.gro" ]; then
    echo "Error: em_cg.gro not found!"
    echo "Please run energy minimization first."
    exit 1
fi

# NVT Warmup Stage 1
echo "Running NVT Warmup Stage 1..."
echo "  - Time step: {self.dt/4.0} ps"
echo "  - Steps: {int(self.nvt_steps * 0.04)}"
{gmx_cmd} grompp -f nvt_warmup1.mdp -c em_cg.gro -p topol.top -o nvt_warmup1.tpr -maxwarn 2
{mdrun_base.format(output_prefix='nvt_warmup1')}
echo "  ✓ NVT Warmup Stage 1 completed"

# NVT Warmup Stage 2
echo "Running NVT Warmup Stage 2..."
echo "  - Time step: {self.dt/2.0} ps"
echo "  - Steps: {int(self.nvt_steps * 0.02)}"
{gmx_cmd} grompp -f nvt_warmup2.mdp -c nvt_warmup1.gro -p topol.top -o nvt_warmup2.tpr -maxwarn 2 -t nvt_warmup1.cpt
{mdrun_base.format(output_prefix='nvt_warmup2')}
echo "  ✓ NVT Warmup Stage 2 completed"

# NVT Production
echo "Running NVT Production..."
echo "  - Time step: {self.dt} ps"
echo "  - Steps: {self.nvt_steps}"
{gmx_cmd} grompp -f nvt_production.mdp -c nvt_warmup2.gro -p topol.top -o nvt_production.tpr -maxwarn 2 -t nvt_warmup2.cpt
{mdrun_base.format(output_prefix='nvt_production')}
echo "  ✓ NVT Production completed"

echo ""
echo "=========================================="
echo "NVT Equilibration completed successfully!"
echo "Final output: nvt_production.gro"
echo "=========================================="
"""
        return script

    def _create_npt_script(self):
        """Create NPT manual run script."""
        gmx_cmd = self.gmx_execute if self.gmx_execute else "gmx"
        
        # Build mdrun command with optional GPU support
        mdrun_base = f"{gmx_cmd} mdrun -v -deffnm {{output_prefix}} -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}"
        if self.use_gpu:
            mdrun_base += " -bonded gpu -update gpu"
        
        # Determine pressure coupling type
        if self.topol == 'slab':
            pressure_info = f"Semi-isotropic pressure coupling (slab topology)"
            pressure_params = f"ref_p = {self.pressure} {self.pressure}"
        else:
            pressure_info = f"Isotropic pressure coupling (grid topology)"
            pressure_params = f"ref_p = {self.pressure}"
        
        script = f"""#!/bin/bash
# Manual NPT equilibration script
# Generated by run_solvent.py
# 
# This script runs NPT equilibration with warmup stages:
# 1. NPT Warmup1: dt={self.dt/4.0} ps, steps={int(self.npt_steps * 0.04)}
# 2. NPT Warmup2: dt={self.dt/2.0} ps, steps={int(self.npt_steps * 0.02)}
# 3. NPT Production: dt={self.dt} ps, steps={self.npt_steps}
# 
# Pressure coupling: {pressure_info}

set -e  # Exit on any error

echo "=========================================="
echo "Starting NPT Equilibration with Warmup"
echo "=========================================="
echo "GROMACS command: {gmx_cmd}"
echo "MPI threads: {self.num_mpi_threads}"
echo "OpenMP threads: {self.num_omp_threads}"
echo "Temperature: {self.temperature} K"
echo "Pressure: {self.pressure} bar"
echo "Pressure coupling: {pressure_info}"
echo "GPU acceleration: {'Enabled' if self.use_gpu else 'Disabled'}"
echo ""

# Check if input file exists
if [ ! -f "nvt_production.gro" ]; then
    echo "Error: nvt_production.gro not found!"
    echo "Please run NVT equilibration first."
    exit 1
fi

# NPT Warmup Stage 1
echo "Running NPT Warmup Stage 1..."
echo "  - Time step: {self.dt/4.0} ps"
echo "  - Steps: {int(self.npt_steps * 0.04)}"
{gmx_cmd} grompp -f npt_warmup1.mdp -c nvt_production.gro -p topol.top -o npt_warmup1.tpr -maxwarn 2 -t nvt_production.cpt
{mdrun_base.format(output_prefix='npt_warmup1')}
echo "  ✓ NPT Warmup Stage 1 completed"

# NPT Warmup Stage 2
echo "Running NPT Warmup Stage 2..."
echo "  - Time step: {self.dt/2.0} ps"
echo "  - Steps: {int(self.npt_steps * 0.02)}"
{gmx_cmd} grompp -f npt_warmup2.mdp -c npt_warmup1.gro -p topol.top -o npt_warmup2.tpr -maxwarn 2 -t npt_warmup1.cpt
{mdrun_base.format(output_prefix='npt_warmup2')}
echo "  ✓ NPT Warmup Stage 2 completed"

# NPT Production
echo "Running NPT Production..."
echo "  - Time step: {self.dt} ps"
echo "  - Steps: {self.npt_steps}"
{gmx_cmd} grompp -f npt_production.mdp -c npt_warmup2.gro -p topol.top -o npt_production.tpr -maxwarn 2 -t npt_warmup2.cpt
{mdrun_base.format(output_prefix='npt_production')}
echo "  ✓ NPT Production completed"

echo ""
echo "=========================================="
echo "NPT Equilibration completed successfully!"
echo "Final output: npt_production.gro"
echo "=========================================="
"""
        return script

    def _create_equilibration_mdps(self):
        """Create MDP files for equilibration with warmup stages."""
        constraints_str = self._get_constraints_string()
        
        # Calculate warmup parameters
        # Warmup 1: 1/4 step size, 4% steps
        dt_warmup1 = self.dt / 4.0
        nvt_steps_warmup1 = int(self.nvt_steps * 0.04)
        npt_steps_warmup1 = int(self.npt_steps * 0.04)
        
        # Warmup 2: 1/2 step size, 2% steps  
        dt_warmup2 = self.dt / 2.0
        nvt_steps_warmup2 = int(self.nvt_steps * 0.02)
        npt_steps_warmup2 = int(self.npt_steps * 0.02)
        
        # Create NVT MDP files (3 stages: warmup1, warmup2, production)
        nvt_warmup1_mdp = self._create_nvt_mdp_template(constraints_str, dt_warmup1, nvt_steps_warmup1)
        nvt_warmup2_mdp = self._create_nvt_mdp_template(constraints_str, dt_warmup2, nvt_steps_warmup2)
        nvt_production_mdp = self._create_nvt_mdp_template(constraints_str, self.dt, self.nvt_steps)
        
        # Create NPT MDP files (3 stages: warmup1, warmup2, production)
        npt_warmup1_mdp = self._create_npt_mdp_template(constraints_str, dt_warmup1, npt_steps_warmup1)
        npt_warmup2_mdp = self._create_npt_mdp_template(constraints_str, dt_warmup2, npt_steps_warmup2)
        npt_production_mdp = self._create_npt_mdp_template(constraints_str, self.dt, self.npt_steps)
        
        # Write all MDP files
        self._write_mdp('nvt_warmup1.mdp', nvt_warmup1_mdp, self.solvent_dir)
        self._write_mdp('nvt_warmup2.mdp', nvt_warmup2_mdp, self.solvent_dir)
        self._write_mdp('nvt_production.mdp', nvt_production_mdp, self.solvent_dir)
        self._write_mdp('npt_warmup1.mdp', npt_warmup1_mdp, self.solvent_dir)
        self._write_mdp('npt_warmup2.mdp', npt_warmup2_mdp, self.solvent_dir)
        self._write_mdp('npt_production.mdp', npt_production_mdp, self.solvent_dir)
        
        # Also create legacy files for backward compatibility
        self._write_mdp('nvt.mdp', nvt_production_mdp, self.solvent_dir)
        self._write_mdp('npt.mdp', npt_production_mdp, self.solvent_dir)

    def _create_nvt_mdp_template(self, constraints_str, dt, nsteps):
        """Create NVT MDP template with specified parameters."""
        return f"""constraints         = {constraints_str}
integrator          = md
dt                  = {dt}
nsteps              = {nsteps}
nstcomm             = 100
nstxout             = 10000000
nstvout             = 10000000
nstxtcout           = 50000
xtc_grps            = system
nstfout             = 0
nstlog              = 10000
nstenergy           = 10000
ns_type             = grid
nstlist = 20
coulombtype         = reaction-field
rcoulomb = 1.1
epsilon_r           = 15
epsilon_rf = 0
vdwtype             = cut-off
;vdw-modifier = Potential-shift-verlet
vdw-modifier = Force-switch
rvdw-switch = 0.9
rvdw                = 1.1
table-extension     = 1.8

Tcoupl              = v-rescale
tc-grps             = Protein  non-Protein
tau_t               = 1      1
ref_t               = {self.temperature} {self.temperature}
; Energy monitoring
;energygrps          = Protein  non-Protein
; Pressure coupling is off for NVT
Pcoupl              = no
"""

    def _create_npt_mdp_template(self, constraints_str, dt, nsteps):
        """Create NPT MDP template with specified parameters."""
        if self.topol == 'slab':
            # Semi-isotropic pressure coupling for slab topology
            return f"""constraints         = {constraints_str}
integrator          = md
dt                  = {dt}
nsteps              = {nsteps}
nstcomm             = 100
nstxout             = 10000000
nstvout             = 10000000
nstxtcout           = 50000
xtc_grps            = system
nstfout             = 0
nstlog              = 10000
nstenergy           = 10000
ns_type             = grid
nstlist = 20
coulombtype         = reaction-field
rcoulomb = 1.1
epsilon_r           = 15
epsilon_rf = 0
vdwtype             = cut-off
;vdw-modifier = Potential-shift-verlet
vdw-modifier = Force-switch
rvdw-switch = 0.9
rvdw                = 1.1
table-extension     = 1.8

Tcoupl              = v-rescale
tc-grps             = Protein  non-Protein
tau_t               = 1      1
ref_t               = {self.temperature} {self.temperature}
; Energy monitoring
;energygrps          = Protein  non-Protein
; Pressure coupling is on
Pcoupl              = c-rescale
tau_p               = 1
pcoupltype = semiisotropic
compressibility     = 4.5e-5 4.5e-5 
ref_p               = {self.pressure}     {self.pressure}
"""
        else:
            # Isotropic pressure coupling for grid topology (default)
            return f"""constraints         = {constraints_str}
integrator          = md
dt                  = {dt}
nsteps              = {nsteps}
nstcomm             = 100
nstxout             = 10000000
nstvout             = 10000000
nstxtcout           = 50000
xtc_grps            = system
nstfout             = 0
nstlog              = 10000
nstenergy           = 10000
ns_type             = grid
nstlist = 20
coulombtype         = reaction-field
rcoulomb = 1.1
epsilon_r           = 15
epsilon_rf = 0
vdwtype             = cut-off
;vdw-modifier = Potential-shift-verlet
vdw-modifier = Force-switch
rvdw-switch = 0.9
rvdw                = 1.1
table-extension     = 1.8

Tcoupl              = v-rescale
tc-grps             = Protein  non-Protein
tau_t               = 1      1
ref_t               = {self.temperature} {self.temperature}
; Energy monitoring
;energygrps          = Protein  non-Protein
; Pressure coupling is on
Pcoupl              = c-rescale
tau_p               = 1
pcoupltype = isotropic
compressibility     = 4.5e-5
ref_p               = {self.pressure}
"""

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
