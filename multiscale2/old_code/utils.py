import subprocess
import shlex
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command, cwd="."):
    """
    Executes a shell command safely and logs its output.

    Args:
        command (str): The command to execute.
        cwd (str): The working directory to run the command in.

    Returns:
        subprocess.CompletedProcess: The result of the command execution.
    """
    logging.info(f"Executing command in '{cwd}':\n  $ {command}")
    try:
        # Using shlex.split to handle command arguments safely
        result = subprocess.run(
            shlex.split(command),
            cwd=cwd,
            check=True,
            text=True,
            capture_output=True
        )
        if result.stdout:
            logging.info(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            logging.warning(f"STDERR:\n{result.stderr}")
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}")
        logging.error(f"STDOUT:\n{e.stdout}")
        logging.error(f"STDERR:\n{e.stderr}")
        raise