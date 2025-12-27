# Multiscale² (MS²): An Advanced Workflow for Protein Condensate Simulations

Multiscale² (MS²) is an automated, high-throughput simulation pipeline designed to bridge coarse-grained (CG) and all-atom (AA) representations of protein condensates. By integrating CALVADOS, OpenMM, and GROMACS, MS² provides a robust framework for generating, refining, and equilibrating complex biomolecular systems, culminating in production-ready, all-atom models.

## Core Innovations

-   **Automated Workflow Generation**: MS² utilizes a powerful `WorkflowGenerator` to create bespoke Python scripts for each stage of the simulation, tailored to your specific system configuration. This eliminates manual script editing and ensures reproducibility.
-   **Advanced Equilibration Protocol**: The explicit solvent stage features a sophisticated, multi-stage warmup protocol with progressive time step scaling (e.g., `0.25*dt` -> `0.5*dt` -> `dt`). This ensures system stability before production runs.
-   **Intelligent System Adaptation**: The workflow automatically adapts simulation parameters based on your system's topology (`grid` vs. `slab`), intelligently selecting the appropriate pressure coupling (isotropic vs. semi-isotropic).
-   **GPU-Accelerated Equilibration**: NVT and NPT equilibration stages can be offloaded to GPUs (`-bonded gpu -update gpu`), dramatically reducing runtime for these computationally intensive steps.
-   **Fine-grained User Control**: While the initial stages are automated, the final equilibration steps are packaged into executable shell scripts (`run_nvt.sh`, `run_npt.sh`), giving the user full control over launching the most resource-intensive parts of the workflow.

## The MS² Simulation Pipeline

The workflow seamlessly transitions a system from a coarse-grained model to a production-ready, all-atom state in explicit solvent.

1.  **`Stage 1: Coarse-Grained Dynamics (CALVADOS)`**
    *   Initial simulation of condensate formation using the CALVADOS CG model.

2.  **`Stage 2: Atomistic Reconstruction (Backmapping)`**
    *   Conversion of the final CG structure into an initial all-atom representation.

3.  **`Stage 3: Steric Refinement (OpenMM)`**
    *   Energy minimization using OpenMM to resolve steric clashes and refine the backmapped structure.

4.  **`Stage 4: Explicit Solvation & Equilibration (GROMACS)`**
    *   Solvation of the refined structure, followed by a multi-stage NVT/NPT equilibration protocol to prepare the system for production MD.

## Quick Start

### Step 0: Generate a New Workflow
For a new system, first generate the necessary workflow scripts.

```bash
# Example: Generate scripts for a new project
python -c "
from multiscale2.workflow_generator import WorkflowGenerator
generator = WorkflowGenerator('path/to/your/new_config.yaml')
generator.generate_all_scripts('path/to/your/new_workflow_dir')
"
```

### Step 1: Configure the System
Define your system and simulation parameters in the `config.yaml` within your workflow directory.

```yaml
# In your_workflow_dir/config.yaml
cg_calvados:
  topol: "slab"  # grid or slab

gromacs_explicit_solvent_build_and_equilibration:
  dt: 0.002
  NVT_steps: 500000
  NPT_steps: 500000
  use_gpu: True  # Set to True to enable GPU offloading
```

### Step 2: Execute Automated Stages
Run the auto-generated Python scripts in sequence to process the system up to the equilibration stage.

```bash
# Navigate to your workflow directory
cd path/to/your/new_workflow_dir

python run_cg.py       # Stage 1
python run_backmap.py  # Stage 2
python run_openmm.py   # Stage 3
python run_solvent.py  # Stage 4 (prepares for equilibration)
```

### Step 3: Launch Manual Equilibration
The final step is to run the generated shell scripts to perform the NVT and NPT equilibration.

```bash
# Navigate to the explicit solvent output directory
cd output_solvent/explicit_solvent/

# Execute the equilibration workflows
./run_nvt.sh
./run_npt.sh
```


## Requirements

-   Python 3.7+
-   GROMACS 2020+
-   OpenMM 7.6+
-   CALVADOS
-   CUDA Toolkit (for GPU acceleration)
-   CG2ALL