from .utils import run_command

class Backmapper:
    """
    A wrapper for backmapping tools like cg2all.
    """
    def __init__(self, config):
        self.config = config
        self.tool = config.get('tool', 'cg2all')
        self.cg_bead = config.get('cg_bead_name', 'CA')
        self.fix = "--fix" if config.get('fix_sidechains', True) else ""

    def reconstruct(self, input_cg_pdb, output_aa_pdb):
        """
        Performs the reconstruction from CG to AA.
        """
        if self.tool == 'cg2all':
            # This command is based on 'prepare_for_backmap.py'
            # Note: A real implementation might need more sophisticated setup
            # like centering the box first.
            command = f"cg2all -p {input_cg_pdb} -o {output_aa_pdb} --cg {self.cg_bead} {self.fix}"
            run_command(command)
        else:
            raise NotImplementedError(f"Backmapping tool '{self.tool}' is not supported.")