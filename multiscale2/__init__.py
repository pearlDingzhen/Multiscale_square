#!/usr/bin/env python3
"""
Multiscale2 Package

A multi-stage workflow for simulating protein condensates.

Subpackages:
- src: Coarse-grained simulation framework
- extern: External packages (ms2_calvados, etc.)
"""

from .src import (
    CGSimulationConfig,
    CGComponent,
    ComponentType,
    TopologyType,
    ComputePlatform,
    SimulationParams,
    SimulationResult,
    CGSimulator,
    CalvadosWrapper,
    run_calvados,
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
    
    # CALVADOS wrapper
    'CalvadosWrapper',
    'run_calvados',
]

