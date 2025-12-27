"""
Workflow Stages Package

This package contains templates for each stage of the multiscale simulation workflow.
Each module is a complete, runnable Python script that can be copied and executed independently.

Modules:
- run_calvados.py: Stage 1 - CALVADOS coarse-grained simulation
- run_backmap.py: Stage 2 - Backmapping CG to all-atom
- run_openmm.py: Stage 3 - OpenMM structure optimization
- run_solvent.py: Stage 4 - Explicit solvent preparation

Usage:
    from multiscale2.workflow_stages import run_calvados
    # or copy the file and run directly:
    # python path/to/workflow_stages/run_calvados.py
"""

__version__ = "1.0.0"
__author__ = "Multiscale Square Team"


