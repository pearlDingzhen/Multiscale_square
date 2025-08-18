import os
import shutil
import yaml
from .calvados import sim, cfg

def run_calvados_simulation(config, output_dir, components_path):
    """
    A wrapper to run the CALVADOS simulation as the first stage.
    It prepares a compatible config file from the main config.
    """
    # Extract the 'cg_calvados' section for compatibility
    calvados_config_dict = config.get('cg_calvados', {})
    if not calvados_config_dict:
        raise ValueError("Configuration for 'cg_calvados' not found in config file.")

    # Write a temporary, CALVADOS-compatible config file
    temp_config_path = os.path.join(output_dir, 'cg_config.yaml')
    with open(temp_config_path, 'w') as f:
        yaml.dump(calvados_config_dict, f)

    # Copy the components file to the output directory
    temp_components_path = os.path.join(output_dir, 'components.yaml')
    shutil.copy(components_path, temp_components_path)

    # Run the simulation using the CALVADOS entry point
    sim.run(path=output_dir, fconfig='cg_config.yaml', fcomponents='components.yaml')