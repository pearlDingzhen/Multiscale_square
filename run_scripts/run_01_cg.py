import sys
import os
# Add project root to path to allow importing multiscale2
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from multiscale2.workflow import MultiscaleWorkflow

def main():
    config_path = os.path.join(os.path.dirname(sys.path[0]), '00_input', 'config.yaml')
    workflow = MultiscaleWorkflow(config_path)
    workflow.execute_stage_1_cg()

if __name__ == "__main__":
    main()