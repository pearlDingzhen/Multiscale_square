#!/usr/bin/env python3
"""
Workflow Generator for Multiscale Simulation
Generates run scripts for each stage by copying templates from workflow_stages/ directory.

Architecture:
- Templates are stored as independent Python files in workflow_stages/
- This generator simply copies templates to the target directory
- Each template is a complete, runnable Python script
"""

import os
import shutil
import yaml
from pathlib import Path


class WorkflowGenerator:
    """Generates run scripts by copying templates from workflow_stages/."""
    
    # List of stage scripts to generate (in order)
    STAGES = ['calvados', 'backmap', 'openmm', 'solvent']
    
    def __init__(self, config_path):
        """Initialize with configuration file path."""
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.config_path = config_path
        self.config_dir = os.path.dirname(os.path.abspath(config_path))
        
        # Global task type (IDP or MDP)
        self.task_type = (
            self.config.get('task_type')
            or self.config.get('cg_calvados', {}).get('task_type')
            or 'IDP'
        )
        
        # Find the workflow_stages directory
        self.stages_dir = self._find_stages_dir()
    
    def _find_stages_dir(self) -> Path:
        """Find the workflow_stages directory containing stage templates."""
        import multiscale2
        package_dir = Path(multiscale2.__file__).parent
        stages_dir = package_dir / "workflow_stages"
        
        if stages_dir.exists():
            return stages_dir
        
        # Also check relative to current working directory
        cwd_stages_dir = Path.cwd() / "workflow_stages"
        if cwd_stages_dir.exists():
            return cwd_stages_dir
        
        raise FileNotFoundError(
            f"workflow_stages directory not found.\n"
            f"  - Looked in: {stages_dir}\n"
            f"  - Looked in: {cwd_stages_dir}\n"
            f"  - Please ensure templates exist in the workflow_stages/ directory."
        )
    
    def _get_template_path(self, stage_name: str) -> Path:
        """Get the path to a template file."""
        template_path = self.stages_dir / f"run_{stage_name}.py"
        
        if not template_path.exists():
            raise FileNotFoundError(
                f"Template not found for stage '{stage_name}':\n"
                f"  - Expected: {template_path}\n"
                f"  - Available templates: {list(self.stages_dir.glob('run_*.py'))}"
            )
        
        return template_path
    
    def _copy_template(self, stage_name: str, output_dir: str) -> str:
        """Copy a stage template to the output directory."""
        template_path = self._get_template_path(stage_name)
        output_path = os.path.join(output_dir, f"run_{stage_name}.py")
        
        # Copy the template file
        shutil.copy2(template_path, output_path)
        
        print(f"Generated: {output_path}")
        return output_path
    
    def generate_all_scripts(self, output_dir="."):
        """Generate all run scripts for the workflow."""
        print("=" * 60)
        print("Multiscale Simulation Workflow Generator")
        print("=" * 60)
        print(f"Task type: {self.task_type}")
        print(f"Template directory: {self.stages_dir}")
        print(f"Output directory: {os.path.abspath(output_dir)}")
        print()
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate each stage script
        for stage in self.STAGES:
            self._copy_template(stage, output_dir)
        
        print()
        print("=" * 60)
        print("Workflow scripts generated successfully!")
        print("=" * 60)
        print()
        print("Generated files:")
        for stage in self.STAGES:
            print(f"  - run_{stage}.py")
        print()
        print("Next steps:")
        print("  1. python run_calvados.py")
        print("  2. python run_backmap.py")
        print("  3. python run_openmm.py")
        print("  4. python run_solvent.py")
        print()


def main():
    """Main function to generate workflow scripts."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate multiscale workflow scripts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python workflow_generator.py config.yaml
  python workflow_generator.py config.yaml -o workflow/
        """
    )
    parser.add_argument(
        "config",
        help="Path to configuration YAML file"
    )
    parser.add_argument(
        "--output", "-o",
        default=".",
        help="Output directory for generated scripts (default: current directory)"
    )
    parser.add_argument(
        "--list-templates",
        action="store_true",
        help="List available templates and exit"
    )
    
    args = parser.parse_args()
    
    try:
        generator = WorkflowGenerator(args.config)
        
        if args.list_templates:
            print("Available templates:")
            for stage in generator.STAGES:
                template_path = generator._get_template_path(stage)
                print(f"  - run_{stage}.py -> {template_path}")
            return
        
        generator.generate_all_scripts(args.output)
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
