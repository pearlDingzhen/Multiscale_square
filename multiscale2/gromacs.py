import os
from .utils import run_command
from . import restraint_generator # New module for plumed logic

class GromacsRunner:
    """
    A wrapper for executing common GROMACS workflows.
    """
    def __init__(self, config, working_dir):
        self.config = config
        self.dir = working_dir
        os.makedirs(self.dir, exist_ok=True)

    def _write_mdp(self, filename, content):
        with open(os.path.join(self.dir, filename), 'w') as f:
            f.write(content)

    def _prepare_mdps(self):
        # Based on mdps from prepare_for_backmap.py
        em_steep_mdp = """
integrator  = steep
nsteps      = 5000
emtol       = 1000.0
emstep      = 0.01
"""
        em_cg_mdp = """
integrator  = cg
nsteps      = 5000
emtol       = 100.0
"""
        npt_mdp = f"""
integrator  = md
nsteps      = {self.config['equilibration']['npt_steps']}
dt          = 0.002
coulombtype = PME
vdwtype     = Cut-off
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = {self.config['equilibration']['temperature']} {self.config['equilibration']['temperature']}
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = {self.config['equilibration']['pressure']}
gen_vel     = yes
gen_temp    = {self.config['equilibration']['temperature']}
"""
        md_mdp = f"""
integrator  = md
nsteps      = {self.config['production']['steps']}
dt          = {self.config['production']['dt']}
coulombtype = PME
vdwtype     = Cut-off
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = {self.config['production']['temperature']}
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = {self.config['production']['pressure']}
"""
        self._write_mdp('em_steep.mdp', em_steep_mdp)
        self._write_mdp('em_cg.mdp', em_cg_mdp)
        self._write_mdp('npt.mdp', npt_mdp)
        self._write_mdp('md.mdp', md_mdp)


    def build_and_solvate(self, input_pdb, num_chains):
        """
        Runs pdb2gmx, solvate, and genion.
        Refactors logic from prepare_for_backmap.py
        """
        # This is a highly simplified version. The original script has complex
        # logic involving PCcli and manual topology editing.
        print("Building GROMACS topology (pdb2gmx)...")
        run_command(f"gmx pdb2gmx -f {input_pdb} -o protein.gro -ignh -ff {self.config['forcefield']} -water none", cwd=self.dir)

        print("Editing topology for multiple chains...")
        # Placeholder for manual topology editing logic

        print("Creating box and solvating...")
        run_command(f"gmx editconf -f protein.gro -o protein_box.gro -c -d {self.config['box_distance']} -bt {self.config['box_type']}", cwd=self.dir)
        run_command(f"gmx solvate -cp protein_box.gro -cs -o system_solv.gro -p topol.top", cwd=self.dir)

        print("Adding ions...")
        run_command(f"gmx grompp -f em_steep.mdp -c system_solv.gro -p topol.top -o ions.tpr -maxwarn 2", cwd=self.dir)
        ion_cmd = f"gmx genion -s ions.tpr -o system_solv_ions.gro -p topol.top -pname NA -nname CL -conc {self.config['ion_conc']}"
        if self.config.get('neutralize', True):
            ion_cmd += " -neutral"
        # Use subprocess with input to handle prompts
        subprocess.run(shlex.split(ion_cmd), input='SOL', text=True, cwd=self.dir, check=True)


    def run_equilibration(self, input_gro, input_top, ref_pdb_for_restraint=None):
        """
        Runs minimization and NVT/NPT equilibration.
        """
        self._prepare_mdps()
        shutil.copy(input_gro, os.path.join(self.dir, 'system.gro'))
        shutil.copy(input_top, os.path.join(self.dir, 'topol.top'))

        plumed_cmd = ""
        if self.config.get('equilibration', {}).get('use_contact_restraints'):
            print("Generating PLUMED restraint file...")
            plumed_file = os.path.join(self.dir, 'plumed.dat')
            restraint_generator.generate_plumed_dat(
                ref_pdb=ref_pdb_for_restraint,
                domains=self.config['equilibration']['restraint_domains'],
                num_monomers=self.config['num_chains'], # Assuming num_chains is in config
                output_file=plumed_file
            )
            plumed_cmd = f"-plumed {plumed_file}"

        print("Running energy minimization...")
        run_command(f"gmx grompp -f em_steep.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 2", cwd=self.dir)
        run_command(f"gmx mdrun -v -deffnm em", cwd=self.dir)

        print("Running NPT equilibration...")
        run_command(f"gmx grompp -f npt.mdp -c em.gro -p topol.top -o npt.tpr -maxwarn 2", cwd=self.dir)
        run_command(f"gmx mdrun -v -deffnm npt {plumed_cmd}", cwd=self.dir)

    def run_production(self, input_gro, input_top, input_cpt):
        """
        Runs the final production simulation.
        """
        self._prepare_mdps()
        shutil.copy(input_gro, os.path.join(self.dir, 'system.gro'))
        shutil.copy(input_top, os.path.join(self.dir, 'topol.top'))
        shutil.copy(input_cpt, os.path.join(self.dir, 'state.cpt'))

        print("Running production simulation...")
        run_command(f"gmx grompp -f md.mdp -c system.gro -t state.cpt -p topol.top -o md.tpr -maxwarn 2", cwd=self.dir)
        run_command(f"gmx mdrun -v -deffnm md", cwd=self.dir)