#!/usr/bin/env python3
"""
Multiscale2 CLI

A multi-stage workflow for simulating protein condensates.

    Commands:
        init    Initialize a new configuration template
        cg      Run coarse-grained simulation
        backmap Backmap CG structure to all-atom representation
        info    Display system information
"""

import sys
import click

from .commands import init_command, cg_command, backmap_command, info_command


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    """
    Multiscale2: A multi-stage workflow for simulating protein condensates.
    
    Commands:
        init    Initialize a new configuration template
        cg      Run coarse-grained simulation
        info    Display system and environment information
    
    Available force fields:
        calvados, hps_urry, cocomo, mpipi_recharged
    
    Examples:
        ms2 init my_project
        ms2 init --ff mpipi_recharged --type mdp
        ms2 cg -f config.yaml -ff mpipi_recharged
        ms2 backmap -i TDP43_CG -f config.yaml
        ms2 info
    """
    pass


# Add commands
main.add_command(init_command, 'init')
main.add_command(cg_command, 'cg')
main.add_command(backmap_command, 'backmap')
main.add_command(info_command, 'info')


def cli():
    """Entry point for CLI."""
    sys.exit(main())


if __name__ == '__main__':
    cli()

