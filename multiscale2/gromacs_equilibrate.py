import os
import shutil
import subprocess
import shlex
from pathlib import Path
from .utils import run_command

class IDPGromacsEquilibrator:
    """
    A specialized GROMACS equilibrator for IDP (Intrinsically Disordered Proteins) systems.
    Handles energy minimization and equilibration (NVT/NPT) from OpenMM optimization output.
    """
    
    def __init__(self, config, working_dir):
        """
        Initialize the IDP GROMACS equilibrator.
        
        Args:
            config: Configuration dictionary containing simulation parameters
            working_dir: Working directory for GROMACS simulations
        """
        self.config = config
        self.working_dir = Path(working_dir).resolve()
        self.working_dir.mkdir(parents=True, exist_ok=True)
        
        # Extract key parameters from config
        self.protein_name = config['protein']['name']
        self.nmol = config['protein']['nmol']
        self.task_type = config.get('task_type', 'IDP')
        
        # Extract GROMACS parameters from config
        gromacs_config = config.get('gromacs_explicit_solvent_build_and_equilibration', {})
        self.num_cpus = gromacs_config.get('num_cpus', 1)
        self.num_mpi_threads = gromacs_config.get('num_mpi_threads', 1)
        self.num_omp_threads = gromacs_config.get('num_omp_threads', self.num_cpus)
        self.gmx_execute = gromacs_config.get('gmx_excute', None)  # Note: keeping original typo from config
        self.ions_concentration = gromacs_config.get('ions_concentration', 0.15)
        self.temperature = gromacs_config.get('temperature', 300.0)
        self.pressure = gromacs_config.get('pressure', 1.0)
        self.use_anti_freeze_water = gromacs_config.get('Use_anti_freeze_water', True)
        self.anti_free_water_ratio = gromacs_config.get('Anti_free_water_ratio', 0.1)
        self.nvt_enabled = gromacs_config.get('NVT', True)
        self.nvt_steps = gromacs_config.get('NVT_steps', 1000000)
        self.npt_enabled = gromacs_config.get('NPT', True)
        self.npt_steps = gromacs_config.get('NPT_steps', 1000000)
        
        # Set up subdirectories
        self.em_dir = self.working_dir / "energy_minimization"
        self.eq_dir = self.working_dir / "equilibration"
        
        for dir_path in [self.em_dir, self.eq_dir]:
            dir_path.mkdir(exist_ok=True)

    def _write_mdp(self, filename, content, target_dir=None):
        """Write MDP file content to specified directory."""
        if target_dir is None:
            target_dir = self.em_dir
        mdp_path = target_dir / filename
        with open(mdp_path, 'w') as f:
            f.write(content)
        print(f"✓ Created MDP file: {mdp_path}")

    def get_openmm_outputs(self, openmm_output_dir):
        """
        Step 1: Get conf.gro and PACE.top from OpenMM optimization output.
        
        Args:
            openmm_output_dir: Path to the OpenMM output directory (e.g., output_openmm/openmm)
            
        Returns:
            tuple: (conf_gro_path, pace_top_path) - paths to the input files
        """
        print("="*60)
        print("Step 1: Getting OpenMM Output Files")
        print("="*60)
        
        openmm_dir = Path(openmm_output_dir).resolve()
        
        # Expected file paths
        conf_gro_path = openmm_dir / "conf.gro"
        pace_top_path = openmm_dir / "PACE.top"
        
        # Check if files exist
        if not conf_gro_path.exists():
            raise FileNotFoundError(f"conf.gro not found at: {conf_gro_path}")
        if not pace_top_path.exists():
            raise FileNotFoundError(f"PACE.top not found at: {pace_top_path}")
        
        print(f"✓ Found conf.gro: {conf_gro_path}")
        print(f"✓ Found PACE.top: {pace_top_path}")
        
        # Copy files to energy minimization directory
        em_conf_gro = self.em_dir / "conf.gro"
        em_pace_top = self.em_dir / "PACE.top"
        
        shutil.copy2(conf_gro_path, em_conf_gro)
        shutil.copy2(pace_top_path, em_pace_top)
        
        print(f"✓ Copied files to energy minimization directory:")
        print(f"  - {em_conf_gro}")
        print(f"  - {em_pace_top}")
        
        return em_conf_gro, em_pace_top

    def run_idp_energy_minimization(self, conf_gro_path, pace_top_path):
        """
        Step 2: Run energy minimization for IDP system with explicit solvent.
        
        Args:
            conf_gro_path: Path to the conf.gro file from OpenMM
            pace_top_path: Path to the PACE.top file from OpenMM
            
        Returns:
            tuple: (minimized_gro_path, minimized_top_path) - paths to minimized files
        """
        print("="*60)
        print("Step 2: IDP Energy Minimization with Explicit Solvent")
        print("="*60)
        
        # Create energy minimization MDP files
        self._create_em_mdps()
        
        # Run energy minimization
        print("Running energy minimization...")
        
        # First minimization: steep descent
        print("  - Steep descent minimization...")
        run_command(
            f"gmx grompp -f em_steep.mdp -c conf.gro -p PACE.top -o em_steep.tpr -maxwarn 2",
            cwd=self.em_dir
        )
        run_command(
            f"gmx mdrun -v -deffnm em_steep -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
            cwd=self.em_dir
        )
        
        # Second minimization: conjugate gradient
        print("  - Conjugate gradient minimization...")
        run_command(
            f"gmx grompp -f em_cg.mdp -c em_steep.gro -p PACE.top -o em_cg.tpr -maxwarn 2",
            cwd=self.em_dir
        )
        run_command(
            f"gmx mdrun -v -deffnm em_cg -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
            cwd=self.em_dir
        )
        
        # Final minimized files
        minimized_gro = self.em_dir / "em_cg.gro"
        minimized_top = self.em_dir / "PACE.top"
        
        print(f"✓ Energy minimization completed!")
        print(f"✓ Minimized structure: {minimized_gro}")
        print(f"✓ Topology file: {minimized_top}")
        
        return minimized_gro, minimized_top

    def _create_em_mdps(self):
        """Create MDP files for energy minimization using the provided template."""
        
        # Steep descent minimization
        em_steep_mdp = f"""define              =  
constraints         =  none
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
        
        # Conjugate gradient minimization  
        em_cg_mdp = f"""define              =  
constraints         =  none
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
        
        self._write_mdp('em_steep.mdp', em_steep_mdp, self.em_dir)
        self._write_mdp('em_cg.mdp', em_cg_mdp, self.em_dir)

    def _create_equilibration_mdps(self):
        """Create MDP files for NVT and NPT equilibration."""
        
        # NVT equilibration
        nvt_mdp = f"""define              =  
constraints         =  h-bonds
integrator          =  md
nsteps              =  {self.nvt_steps}
dt                  =  0.002
lincs_iter          =  8
emtol               =  1000.0
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

tcoupl              =  V-rescale
tc-grps             =  Protein Non-Protein
tau_t               =  0.1 0.1
ref_t               =  {self.temperature} {self.temperature}
pcoupl              =  no
gen_vel             =  yes
gen_temp            =  {self.temperature}

nstcomm             =  100
table-extension     =  10
"""
        
        # NPT equilibration
        npt_mdp = f"""define              =  
constraints         =  h-bonds
integrator          =  md
nsteps              =  {self.npt_steps}
dt                  =  0.002
lincs_iter          =  8
emtol               =  1000.0
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

tcoupl              =  V-rescale
tc-grps             =  Protein Non-Protein
tau_t               =  0.1 0.1
ref_t               =  {self.temperature} {self.temperature}
pcoupl              =  Parrinello-Rahman
pcoupltype          =  isotropic
tau_p               =  2.0
ref_p               =  {self.pressure}
compressibility     =  4.5e-5
gen_vel             =  no

nstcomm             =  100
table-extension     =  10
"""
        
        self._write_mdp('nvt.mdp', nvt_mdp, self.eq_dir)
        self._write_mdp('npt.mdp', npt_mdp, self.eq_dir)

    def run_equilibration(self, minimized_gro_path, minimized_top_path):
        """
        Step 3: Run NVT and NPT equilibration.
        
        Args:
            minimized_gro_path: Path to the minimized structure from energy minimization
            minimized_top_path: Path to the topology file
            
        Returns:
            tuple: (equilibrated_gro_path, equilibrated_top_path) - paths to equilibrated files
        """
        print("="*60)
        print("Step 3: NVT and NPT Equilibration")
        print("="*60)
        
        # Create equilibration MDP files
        self._create_equilibration_mdps()
        
        # Copy minimized files to equilibration directory
        eq_gro = self.eq_dir / "system.gro"
        eq_top = self.eq_dir / "topol.top"
        
        shutil.copy2(minimized_gro_path, eq_gro)
        shutil.copy2(minimized_top_path, eq_top)
        
        print(f"✓ Copied minimized files to equilibration directory:")
        print(f"  - {eq_gro}")
        print(f"  - {eq_top}")
        
        # Run NVT equilibration if enabled
        if self.nvt_enabled:
            print(f"\nRunning NVT equilibration ({self.nvt_steps} steps)...")
            run_command(
                f"gmx grompp -f nvt.mdp -c system.gro -p topol.top -o nvt.tpr -maxwarn 2",
                cwd=self.eq_dir
            )
            run_command(
                f"gmx mdrun -v -deffnm nvt -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
                cwd=self.eq_dir
            )
            print("✓ NVT equilibration completed!")
            
            # Use NVT output for NPT
            nvt_gro = self.eq_dir / "nvt.gro"
            nvt_cpt = self.eq_dir / "nvt.cpt"
        else:
            print("NVT equilibration disabled, using minimized structure for NPT")
            nvt_gro = eq_gro
            nvt_cpt = None
        
        # Run NPT equilibration if enabled
        if self.npt_enabled:
            print(f"\nRunning NPT equilibration ({self.npt_steps} steps)...")
            npt_cmd = f"gmx grompp -f npt.mdp -c {nvt_gro.name} -p topol.top -o npt.tpr -maxwarn 2"
            if nvt_cpt and nvt_cpt.exists():
                npt_cmd += f" -t {nvt_cpt.name}"
            
            run_command(npt_cmd, cwd=self.eq_dir)
            run_command(
                f"gmx mdrun -v -deffnm npt -ntmpi {self.num_mpi_threads} -ntomp {self.num_omp_threads}",
                cwd=self.eq_dir
            )
            print("✓ NPT equilibration completed!")
            
            # Final equilibrated files
            equilibrated_gro = self.eq_dir / "npt.gro"
            equilibrated_cpt = self.eq_dir / "npt.cpt"
        else:
            print("NPT equilibration disabled")
            equilibrated_gro = nvt_gro
            equilibrated_cpt = nvt_cpt
        
        print(f"✓ Equilibration completed!")
        print(f"✓ Final structure: {equilibrated_gro}")
        print(f"✓ Final topology: {eq_top}")
        if equilibrated_cpt:
            print(f"✓ Checkpoint file: {equilibrated_cpt}")
        
        return equilibrated_gro, eq_top, equilibrated_cpt

    def run_idp_equilibration_workflow(self, openmm_output_dir):
        """
        Main method to run the complete IDP GROMACS equilibration workflow.
        
        Args:
            openmm_output_dir: Path to the OpenMM output directory
            
        Returns:
            dict: Paths to final output files
        """
        print("="*80)
        print(f"Starting IDP GROMACS Equilibration Workflow for {self.protein_name}")
        print(f"Number of molecules: {self.nmol}")
        print(f"Task type: {self.task_type}")
        print(f"Temperature: {self.temperature} K")
        print(f"Pressure: {self.pressure} bar")
        print(f"CPUs: {self.num_cpus}")
        print("="*80)
        
        # Step 1: Get OpenMM outputs
        conf_gro, pace_top = self.get_openmm_outputs(openmm_output_dir)
        
        # Step 2: Energy minimization
        minimized_gro, minimized_top = self.run_idp_energy_minimization(conf_gro, pace_top)
        
        # Step 3: Equilibration (NVT/NPT)
        equilibrated_gro, equilibrated_top, equilibrated_cpt = self.run_equilibration(minimized_gro, minimized_top)
        
        print("="*80)
        print("IDP GROMACS Equilibration Workflow Complete!")
        print("="*80)
        print(f"✓ Minimized structure: {minimized_gro}")
        print(f"✓ Equilibrated structure: {equilibrated_gro}")
        print(f"✓ Topology file: {equilibrated_top}")
        if equilibrated_cpt:
            print(f"✓ Checkpoint file: {equilibrated_cpt}")
        print()
        print("System is now ready for production simulation!")
        print("="*80)
        
        return {
            'minimized_gro': minimized_gro,
            'equilibrated_gro': equilibrated_gro,
            'equilibrated_top': equilibrated_top,
            'equilibrated_cpt': equilibrated_cpt,
            'working_dir': self.working_dir
        }