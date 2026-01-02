#!/usr/bin/env python3
"""
Coarse-Grained Simulation Module

A unified interface for coarse-grained molecular dynamics simulations.

Available runners:
- run_calvados: CALVADOS force field
- run_hps: HPS-Urry force field (OpenABC)
- run_cocomo: COCOMO force field
- run_openmpipi: OpenMpipi force field

Usage:
    from multiscale2.src import CGSimulationConfig, CGSimulator
    
    config = CGSimulationConfig.from_yaml("config.yaml")
    sim = CGSimulator(config)
    sim.setup("output/")
    sim.run_calvados(gpu_id=0)
"""

from .cg import (
    CGSimulationConfig,
    CGComponent,
    ComponentType,
    TopologyType,
    ComputePlatform,
    SimulationParams,
    SimulationResult,
    CGSimulator,
    ChainLabel,
)

from .calvados_wrapper import (
    CalvadosWrapper,
    run_calvados,
)

from .cocomo2_creator import (
    COCOMO,
)

from .pdb_tool import (
    ChainLabel,
    extract_coordinates_from_pdb,
)

__all__ = [
    # Configuration
    'CGSimulationConfig',
    'CGComponent',
    'ComponentType',
    'TopologyType',
    'ComputePlatform',
    'SimulationParams',
    'SimulationResult',

    # Simulator
    'CGSimulator',

    # PDB Tools
    'ChainLabel',
    'extract_coordinates_from_pdb',

    # CALVADOS wrapper
    'CalvadosWrapper',
    'run_calvados',

    # COCOMO creator
    'COCOMO',
]
