#!/usr/bin/env python3
"""
CLI Commands Module

Implements the command-line interface for multiscale2.
"""

import os
import sys
from pathlib import Path
from typing import Optional, List

import click

from ..src import CGSimulationConfig, CGComponent, ComponentType, TopologyType


# Available force fields
FORCE_FIELDS = ['calvados', 'hps_urry', 'cocomo', 'mpipi_recharged']

# Mapping from CLI force field names to internal runner method names
FORCE_FIELD_TO_RUNNER = {
    'calvados': 'calvados',
    'hps_urry': 'hps',  # CLI uses hps_urry, but method is run_hps
    'cocomo': 'cocomo',
    'mpipi_recharged': 'mpipi_recharged',
}


def validate_force_field(ctx, param, value):
    """Validate force field name."""
    if value and value.lower() not in FORCE_FIELDS:
        raise click.BadParameter(
            f"Unknown force field '{value}'. Available: {', '.join(FORCE_FIELDS)}"
        )
    return value.lower() if value else 'calvados'


# =============================================================================
# Init Command
# =============================================================================

@click.command('init', context_settings={'help_option_names': ['-h', '--help']})
@click.argument('name', type=str, required=False)
@click.option(
    '--ff', '-ff',
    type=click.Choice(FORCE_FIELDS),
    default='calvados',
    help=f'Force field. Available: {", ".join(FORCE_FIELDS)} (default: calvados)'
)
@click.option(
    '--type', '-t',
    type=click.Choice(['idp', 'mdp', 'mixed']),
    default='idp',
    help='Component type (default: idp)'
)
@click.option(
    '--topol', '-tp',
    type=click.Choice(['cubic', 'slab']),
    default='cubic',
    help='Topology type (default: cubic)'
)
@click.option(
    '--output', '-o',
    type=click.Path(),
    default='.',
    help='Output directory (default: current directory)'
)
@click.option(
    '--nmol', '-n',
    type=int,
    default=10,
    help='Number of molecules per component (default: 10)'
)
def init_command(name: str, ff: str, type: str, topol: str, output: str, nmol: int):
    """
    Initialize a new CG simulation configuration template.
    
    NAME is the system name (optional, defaults to 'my_simulation').
    
    Examples:
        ms2 init                    # Creates my_simulation.yaml
        ms2 init my_project         # Creates my_project.yaml
        ms2 init --type mdp         # MDP template
        ms2 init --topol slab       # SLAB topology template
    """
    # Set default name if not provided
    if not name:
        name = 'my_simulation'
    
    # Build component list based on type
    components = []
    
    if type == 'idp':
        components.append({
            'name': 'protein_A',
            'type': 'IDP',
            'nmol': nmol,
            'ffasta': 'input/protein_A.fasta',
        })
        component_note = "IDP (requires FASTA file)"
    elif type == 'mdp':
        components.append({
            'name': 'protein_A',
            'type': 'MDP',
            'nmol': nmol,
            'fpdb': 'input/protein_A.pdb',
            'fdomains': 'input/protein_A_domains.yaml',
            'restraint': True,
        })
        component_note = "MDP (requires PDB and domains files)"
    else:  # mixed
        components.extend([
            {
                'name': 'protein_A',
                'type': 'IDP',
                'nmol': nmol,
                'ffasta': 'input/protein_A.fasta',
            },
            {
                'name': 'protein_B',
                'type': 'MDP',
                'nmol': nmol,
                'fpdb': 'input/protein_B.pdb',
                'fdomains': 'input/protein_B_domains.yaml',
                'restraint': True,
            },
        ])
        component_note = "Mixed IDP + MDP"
    
    # Create config
    config = CGSimulationConfig(
        system_name=name,
        box=[25.0, 25.0, 30.0],
        temperature=310.0,
        ionic=0.15,
        topol=TopologyType(topol),
        components=[CGComponent.from_dict(c) for c in components],
    )
    
    # Output path
    output_path = Path(output)
    config_file = output_path / f'{name}.yaml'
    
    # Ensure output directory exists
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save config
    config.to_yaml(str(config_file))
    
    click.echo(f"\n{'=' * 60}")
    click.echo(f"Configuration template created successfully!")
    click.echo(f"{'=' * 60}")
    click.echo(f"\n  File: {config_file}")
    click.echo(f"  System: {name}")
    click.echo(f"  Force field: {ff}")
    click.echo(f"  Topology: {topol}")
    click.echo(f"  Components: {len(components)} ({component_note})")
    click.echo(f"  Molecules per component: {nmol}")
    
    click.echo(f"\n  Next steps:")
    click.echo(f"    1. Edit {config_file}")
    click.echo(f"    2. Add your input files (FASTA/PDB) to 'input/' directory")
    click.echo(f"    3. Run: ms2 cg -f {config_file} -ff {ff}")
    click.echo()


# =============================================================================
# CG Command
# =============================================================================

@click.command('cg', context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--input-file', '-f',
    type=click.Path(exists=True),
    help='Input configuration file (YAML)',
)
@click.option(
    '--force-field', '-ff',
    type=str,
    default='calvados',
    callback=validate_force_field,
    help=f'Force field to use. Available: {", ".join(FORCE_FIELDS)} (default: calvados)',
)
@click.option(
    '--output-dir', '-o',
    type=click.Path(),
    help='Output directory (default: same as input file directory)',
)
@click.option(
    '--dry-run', '-d',
    is_flag=True,
    default=False,
    help='Generate config files without running simulation',
)
@click.option(
    '--overwrite', '-w',
    is_flag=True,
    default=False,
    help='Overwrite existing output directory',
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    default=False,
    help='Enable verbose output',
)
@click.option(
    '--gpu-id', '-g',
    type=int,
    default=0,
    help='GPU device ID (default: 0)',
)
@click.option(
    '--continue-from', '-c',
    type=click.Path(exists=True),
    default=None,
    help='Continue simulation from a PDB coordinate file (Calvados only)',
)
def cg_command(
    input_file: str,
    force_field: str,
    output_dir: Optional[str],
    dry_run: bool,
    overwrite: bool,
    verbose: bool,
    gpu_id: int,
    continue_from: str,
):
    """
    Run coarse-grained simulation.
    
    INPUT_FILE is the configuration YAML file.
    
    Examples:
        ms2 cg -f config.yaml                 # Run with calvados (GPU 0)
        ms2 cg -f config.yaml -g 1            # Run on GPU 1
        ms2 cg -f config.yaml -ff mpipi_recharged  # Run with Mpipi-Recharged
        ms2 cg -f config.yaml -ff cocomo      # Run with COCOMO
        ms2 cg -f config.yaml --dry-run       # Generate config only
        ms2 cg -f config.yaml -o ./results    # Custom output directory
        ms2 cg -f config.yaml -c coords.pdb   # Continue from PDB (Calvados only)
    """
    # Handle input file
    if input_file is None:
        click.echo(click.style("Error: No input file specified.", fg='red'))
        click.echo("Use 'ms2 cg -f config.yaml' or 'ms2 --help' for usage.")
        sys.exit(1)
    
    input_path = Path(input_file).resolve()
    
    # Determine output directory
    if output_dir is None:
        output_dir = str(input_path.parent)
    else:
        output_dir = str(Path(output_dir).resolve())

    click.echo(f"\n{'=' * 60}")
    click.echo(f"CG Simulation")
    click.echo(f"{'=' * 60}")

    click.echo(f"\n  Input file: {input_path}")
    click.echo(f"  Force field: {force_field}")
    click.echo(f"  GPU ID: {gpu_id}")
    click.echo(f"  Mode: {'Dry run' if dry_run else 'Full simulation'}")
    if continue_from:
        click.echo(f"  Continue from: {continue_from}")

    # Load configuration first to get system_name
    click.echo(f"\n[1/4] Loading configuration...")
    try:
        config = CGSimulationConfig.from_yaml(str(input_path))
        click.echo(f"  ✓ Loaded: {config.system_name}")
        click.echo(f"  ✓ Components: {len(config.components)}")
        click.echo(f"  ✓ Total molecules: {config.total_molecules()}")
    except Exception as e:
        click.echo(f"  ✗ Failed to load configuration: {e}")
        sys.exit(1)

    # Calculate the final output directory with _CG suffix
    final_output_dir = os.path.join(output_dir, f"{config.system_name}_CG")
    click.echo(f"  Output dir: {final_output_dir}")
    
    # Validate configuration
    click.echo(f"\n[2/4] Validating configuration...")
    errors = config.validate()
    if errors:
        click.echo(f"  ✗ Validation failed:")
        for error in errors:
            click.echo(f"    - {error}")
        sys.exit(1)
    click.echo(f"  ✓ Configuration valid")
    
    # Validate --continue-from flag
    if continue_from and force_field != 'calvados':
        click.echo(f"  ✗ --continue-from is only supported with Calvados force field.")
        click.echo(f"    Current force field: {force_field}")
        sys.exit(1)

    # Setup simulator with final output directory
    click.echo(f"\n[3/4] Setting up simulation...")
    try:
        from ..src import CGSimulator

        sim = CGSimulator(config)

        # Setup output directory (with _CG suffix)
        sim.setup(final_output_dir, overwrite=overwrite)
        click.echo(f"  ✓ Output directory: {final_output_dir}")

    except FileExistsError:
        click.echo(f"  ✗ Output directory exists: {final_output_dir}")
        click.echo(f"    Use --overwrite to replace.")
        sys.exit(1)
    except Exception as e:
        click.echo(f"  ✗ Setup failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Run simulation or just generate config
    if dry_run:
        click.echo(f"\n[4/4] Dry run - generating configuration files...")

        if force_field == 'calvados':
            # 使用 CGSimulator 的目录逻辑生成配置
            try:
                from ..src.calvados_wrapper import CalvadosWrapper

                # 准备 raw 目录（与 run_calvados 相同的目录结构）
                raw_dir = os.path.join(final_output_dir, 'raw')

                # 写入配置文件到 raw 目录
                wrapper = CalvadosWrapper(config)
                wrapper._write_to_dir(raw_dir)

                click.echo(f"  ✓ Configuration files generated")
                click.echo(f"  Output: {final_output_dir}")
                click.echo(f"  - config.yaml")
                click.echo(f"  - components.yaml")
                click.echo(f"\n  To run simulation, remove --dry-run flag.")
            except Exception as e:
                click.echo(f"  ✗ Failed to generate config: {e}")
                if verbose:
                    import traceback
                    traceback.print_exc()
                sys.exit(1)
        else:
            # 其他力场：先说明需要用 CALVADOS 生成结构
            click.echo(f"  Force field: {force_field}")
            click.echo(f"\n  For non-CALVADOS force fields:")
            click.echo(f"    1. CALVADOS will first generate initial structure")
            click.echo(f"    2. Then the selected force field will run simulation")
            click.echo(f"\n  To run full simulation, remove --dry-run flag.")
    else:
        click.echo(f"\n[4/4] Running {force_field.upper()} simulation...")

        # Call appropriate runner based on force field
        # Map CLI force field name to internal runner method name
        runner_name = FORCE_FIELD_TO_RUNNER.get(force_field, force_field)
        runner_method = f'run_{runner_name}'
        if not hasattr(sim, runner_method):
            click.echo(f"  ✗ Force field '{force_field}' not yet implemented.")
            sys.exit(1)

        try:
            # Call runner method with gpu_id and continue_from parameters
            # For mpipi_recharged, defaults to use_gmx_insert=True, gmx_radius=0.35
            # continue_from is only supported for calvados
            if force_field == 'calvados':
                result = getattr(sim, runner_method)(gpu_id=gpu_id, continue_from=continue_from)
            else:
                result = getattr(sim, runner_method)(gpu_id=gpu_id)

            # Use the result output directory
            actual_output = result.output_dir if result.output_dir else final_output_dir

            if result.success:
                click.echo(f"  ✓ Simulation completed")
                click.echo(f"  Output: {actual_output}")

                if result.trajectory and os.path.exists(result.trajectory):
                    size_mb = os.path.getsize(result.trajectory) / 1024 / 1024
                    click.echo(f"  Trajectory: {os.path.basename(result.trajectory)} ({size_mb:.1f} MB)")
            else:
                click.echo(f"  ✗ Simulation failed:")
                for error in result.errors:
                    click.echo(f"    - {error}")
                sys.exit(1)

        except Exception as e:
            click.echo(f"  ✗ Simulation error: {e}")
            if verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
    
    click.echo()


# =============================================================================
# Backmap Command
# =============================================================================

@click.command('backmap', context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--input', '-i',
    type=click.Path(exists=True),
    required=True,
    help='Input: CG output directory (ms2 cg) or PDB file (user provided)',
)
@click.option(
    '--input-file', '-f',
    type=click.Path(exists=True),
    help='Configuration YAML (same as CG stage, same flag as ms2 cg)',
)
@click.option(
    '--output', '-o',
    type=click.Path(),
    help='Output directory (default: {system_name}_backmap)',
)
@click.option(
    '--device', '-d',
    type=click.Choice(['cpu', 'cuda']),
    help='Device for backmapping (overrides config, default: cpu)',
)
@click.option(
    '--model-type', '-m',
    type=click.Choice(['ResidueBasedModel', 'CalphaBasedModel', 'auto']),
    help='CG model type (overrides config, default: auto-detect)',
)
def backmap_command(input, input_file, output, device, model_type):
    """Backmap CG structure to all-atom representation.
    
    Uses the same config.yaml as CG stage, with optional backmap section.
    
    Examples:
        # ms2 cg output (with explicit config.yaml)
        ms2 backmap -i TDP43_CG -f config.yaml
        
        # ms2 cg output (auto-find config.yaml in pwd)
        ms2 backmap -i TDP43_CG
        
        # user provided PDB (requires config.yaml)
        ms2 backmap -i my_structure.pdb -f config.yaml
        
        # Override device and output
        ms2 backmap -i TDP43_CG -f config.yaml -d cuda -o ./results
    """
    from ..src.backmap import BackmapSimulator, BackmapConfig
    from ..src import CGSimulationConfig
    
    click.echo(f"\n{'=' * 60}")
    click.echo(f"Backmap CG to All-Atom")
    click.echo(f"{'=' * 60}")
    
    click.echo(f"\n  Input: {input}")
    if input_file:
        click.echo(f"  Config: {input_file}")
    if output:
        click.echo(f"  Output: {output}")
    if device:
        click.echo(f"  Device: {device}")
    if model_type:
        click.echo(f"  Model: {model_type}")
    
    # 创建 backmap config（CLI 参数优先）
    backmap_config = BackmapConfig()
    if device:
        backmap_config.device = device
    if model_type and model_type != 'auto':
        backmap_config.model_type = model_type
    if output:
        backmap_config.output_dir = output
    
    # 创建 simulator
    simulator = BackmapSimulator(backmap_config=backmap_config)
    
    # 执行 backmap（output_dir 参数会覆盖 backmap_config.output_dir）
    click.echo(f"\n[1/2] Preparing input...")
    try:
        result = simulator.run(input_path=input, config_path=input_file, output_dir=output)
        
        if result.success:
            click.echo(f"  ✓ Backmap completed")
            click.echo(f"  Input PDB: {result.input_pdb}")
            click.echo(f"  Output PDB: {result.output_pdb}")
            click.echo(f"  Model type: {result.model_type}")
            click.echo(f"\n  ✓ Success!")
        else:
            click.echo(f"  ✗ Backmap failed:")
            for error in result.errors:
                click.echo(f"    - {error}")
            sys.exit(1)
            
    except Exception as e:
        click.echo(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    click.echo()


# =============================================================================
# Pace-Opt Command
# =============================================================================

@click.command('pace-opt', context_settings={'help_option_names': ['-h', '--help']})
@click.option(
    '--input', '-i',
    type=click.Path(exists=True),
    required=True,
    help='Input: backmap output directory (ms2 backmap) or PDB file (user provided)',
)
@click.option(
    '--input-file', '-f',
    type=click.Path(exists=True),
    help='Configuration YAML (same as CG stage, same flag as ms2 cg)',
)
@click.option(
    '--output', '-o',
    type=click.Path(),
    help='Output directory (default: {system_name}_pace_opt)',
)
@click.option(
    '--skip-gromacs',
    is_flag=True,
    default=False,
    help='Skip GROMACS optimization step',
)
@click.option(
    '--device', '-d',
    type=click.Choice(['cpu', 'cuda', 'opencl']),
    default='cuda',
    help='Device for OpenMM (default: cuda)',
)
@click.option(
    '--box-resize',
    type=str,
    help='Box resize dimensions as comma-separated values (e.g., "23.0,17.0,47.0")',
)
@click.option(
    '-l', '--level',
    type=click.Choice(['high', 'medium', 'low']),
    default='medium',
    help='Optimization level: high (7 steps), medium (5 steps), low (3 steps). Default: medium',
)
def pace_opt_command(input, input_file, output, skip_gromacs, device, box_resize, level):
    """Optimize backmapped structure with PACE force field.

    Uses multi-step optimization: Gaussian repulsion → Softcore (multiple steps) → Standard force field.

    Optimization levels:
        high:   7 steps (0.2 → 0.15 → 0.1 → 0.075 → 0.05 → 0.025 → 0.01)
        medium: 5 steps (0.2 → 0.1 → 0.05 → 0.025 → 0.01) - default
        low:    3 steps (0.15 → 0.05 → 0.025)

    Examples:
        # From ms2 backmap output
        ms2 pace-opt -i TDP43_backmap -f config.yaml

        # From user provided PDB
        ms2 pace-opt -i my_structure.pdb -f config.yaml

        # With box resize
        ms2 pace-opt -i TDP43_backmap -f config.yaml --box-resize "23.0,17.0,47.0"

        # Use faster optimization (medium mode)
        ms2 pace-opt -i TDP43_backmap -f config.yaml --optimization-mode medium
        
        # Skip GROMACS optimization
        ms2 pace-opt -i TDP43_backmap -f config.yaml --skip-gromacs
        
        # Use CPU instead of CUDA
        ms2 pace-opt -i TDP43_backmap -f config.yaml -d cpu
    """
    from ..src.pace_opt import PaceOptSimulator, PaceOptConfig
    from ..src import CGSimulationConfig
    
    click.echo(f"\n{'=' * 60}")
    click.echo(f"PACE Force Field Optimization")
    click.echo(f"{'=' * 60}")
    
    click.echo(f"\n  Input: {input}")
    if input_file:
        click.echo(f"  Config: {input_file}")
    if output:
        click.echo(f"  Output: {output}")
    click.echo(f"  Device: {device.upper()}")
    if box_resize:
        click.echo(f"  Box resize: {box_resize}")
    if skip_gromacs:
        click.echo(f"  GROMACS optimization: Skipped")
    
    # Load configuration
    config = None
    if input_file:
        config = CGSimulationConfig.from_yaml(input_file)
    
    # Parse box resize if provided
    box_resize_dims = None
    if box_resize:
        try:
            box_resize_dims = [float(x.strip()) for x in box_resize.split(',')]
            if len(box_resize_dims) != 3:
                raise ValueError("Box resize must have 3 dimensions")
        except ValueError as e:
            click.echo(f"  ✗ Invalid box resize format: {e}")
            sys.exit(1)
    
    # Create pace_opt_config
    pace_opt_config = PaceOptConfig()
    pace_opt_config.box_resize_enabled = box_resize_dims is not None
    pace_opt_config.box_resize_dimensions = box_resize_dims
    pace_opt_config.gromacs_enabled = not skip_gromacs
    pace_opt_config.platform = device.upper()
    pace_opt_config.set_optimization_mode(level)
    
    # Create simulator
    simulator = PaceOptSimulator(config=config, pace_opt_config=pace_opt_config)
    
    # Run optimization
    try:
        result = simulator.run(input_path=input, config_path=input_file, output_dir=output)
        
        if result.success:
            click.echo(f"\n  ✓ PACE optimization completed")
            click.echo(f"  Input PDB: {result.input_pdb}")
            click.echo(f"  Output PDB: {result.output_pdb}")
            if result.intermediate_files:
                click.echo(f"  Intermediate files:")
                for f in result.intermediate_files:
                    click.echo(f"    - {f}")
            click.echo(f"\n  ✓ Success!")
        else:
            click.echo(f"  ✗ PACE optimization failed:")
            for error in result.errors:
                click.echo(f"    - {error}")
            sys.exit(1)
            
    except Exception as e:
        click.echo(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    click.echo()


# =============================================================================
# Info Command
# =============================================================================

@click.command('info', context_settings={'help_option_names': ['-h', '--help']})
def info_command():
    """
    Display system and environment information.
    """
    click.echo(f"\n{'=' * 60}")
    click.echo(f"Multiscale2 Environment")
    click.echo(f"{'=' * 60}")
    
    # Python version
    click.echo(f"\n  Python: {sys.version}")
    
    # Available force fields
    click.echo(f"\n  Available force fields:")
    for ff in FORCE_FIELDS:
        click.echo(f"    - {ff}")
    
    # Check GPU availability
    try:
        import torch
        cuda_available = torch.cuda.is_available()
        if cuda_available:
            click.echo(f"\n  CUDA: Available")
            click.echo(f"    Devices: {torch.cuda.device_count()}")
        else:
            click.echo(f"\n  CUDA: Not available")
    except ImportError:
        click.echo(f"\n  CUDA: Unknown (torch not installed)")
    
    click.echo()

