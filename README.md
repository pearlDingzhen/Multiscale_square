# Multiscale2 (MS2) Simulation Workflow

This project provides a comprehensive, multi-stage workflow for simulating protein condensates, progressing from coarse-grained (CG) models to full all-atom (AA) representations. It uses CALVADOS for the initial CG simulation and GROMACS for subsequent PACE and AA simulations.

## Project Structure

- **/multiscale2**: The core Python package containing all the logic for the workflow.
  - **/calvados**: The original CALVADOS source code for CG simulations.
  - **workflow.py**: The main orchestrator that drives the entire process.
  - **gromacs.py**: A powerful wrapper for all GROMACS command-line operations.
  - **backmap.py**: Handles the backmapping from CG to AA.
  - **openmm_refine.py**: Refines structures using OpenMM to resolve clashes.
  - **aa_transition.py**: Manages the transition from PACE to all-atom models.
  - **utils.py**: Helper utilities.

- **/00_input**: Contains all initial input files.
  - **config.yaml**: The master configuration file for the entire workflow. **EDIT THIS FILE** to control your simulation.
  - **components.yaml**: Defines the molecular components (as in CALVADOS).
  - **protein.fasta**: The FASTA sequence of the protein.
  - **reference_structure.pdb**: A reference PDB for generating contact restraints (optional).

- **/run_scripts**: Executable Python scripts to run each stage of the workflow.

- **/01_cg_calvados ... 07_aa_production**: Directories that will be populated with the output of each stage.

## How to Run the Workflow

1.  **Configure**: Modify `00_input/config.yaml` to set up your specific system and simulation parameters.
2.  **Execute Stages Sequentially**: Run the scripts in the `run_scripts` directory in numerical order. Each script corresponds to one stage of the workflow.

    ```bash
    # Example for Stage 1 (assuming you are in the Multiscale2_Project directory)
    python run_scripts/run_01_cg.py

    # Example for Stage 2
    python run_scripts/run_02_backmap.py

    # ... and so on for all stages.
    ```

Each script takes the output of the previous stage as its input, creating a seamless pipeline from coarse-grained to all-atom simulation.# Multiscale_square
