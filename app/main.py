import os
import subprocess
import shutil
from pathlib import Path

import yaml
from dotenv import load_dotenv

from common import get_repo_directory, inspect_csv
from database import DataBase
from constants import INPUT_TABLE
from runtime_variables import (
    EXECUTION_WORKSPACE,
    CONFIG_YAML,
    ODE_APPLICATION_BINARY,
    DATABASE
)

def prepare_execution_dir(directory: str) -> str:
    """Remove from filesystem what is in "directory" and prepare
    the directory for the execution

    Args:
        directory (str): path to the directory where the program will use for
        the temporary and result files in the execution

    Returns:
        str: the absolute path to execution dir
    """
    print("Preparing execution folder...")
    dir_path = Path(directory)
    if dir_path.exists():
        shutil.rmtree(str(dir_path))
    dir_path.mkdir(parents=True)
    out_path = dir_path.joinpath("out")
    out_path.mkdir()
    print(f"Execution folder prepared: \"{directory}\"")
    return str(out_path)

def prepare_config_file(config_filename: str, destination_path: str) -> str:
    """Read the config.yaml for the execution and modify the "outFolder" field
    to manage where the output is going to be saved. The modified config.yaml is saved
    in the "destination_path" directory

    Args:
        config_filename (str): name of the config.yaml for the execution
        destination_path (str): where the new config.yaml will be saved

    Returns:
        str: the absolute path to the new "config.yaml" file
    """
    
    CONFIG_FILE_MODIFIED = "config.yaml"

    print("Preparing configuration file...")
    # Change outFolder name of the config_file
    with open(config_filename, "r") as f:
        yaml_dict: dict = yaml.safe_load(f)
    out_destination_path = Path(destination_path).joinpath("out")
    yaml_dict["outFileService"]["outFolder"] = str(out_destination_path)

    # Write changes to the destination path
    config_tmp_path = Path(destination_path).joinpath(CONFIG_FILE_MODIFIED)
    with open(config_tmp_path, "w") as f:
        yaml.dump(yaml_dict, f)

    print(f"Configuration file prepared: \"{config_tmp_path}\"")
    return str(config_tmp_path)

def register_input_if_needed(config_filename: str, db_connection: DataBase) -> str:
    """Register input if it does not exists in the DataBase and return the input_id

    Args:
        config_filename (str): path to the configuration yaml file
        db_connection (DataBase): DataBase object already initialized

    Returns:
        str: input_id (path where the input csv file is located)
    """
    print("Registering new input in the database...")
    with open(config_filename, "r") as f:
        config_yaml: dict = yaml.safe_load(f)
    id: str = config_yaml["reader"]["filename"]
    mechanism: str = config_yaml["integrator"]["mechanism"]
    mechanism = mechanism.replace(".yaml","")
    nsp, systems = inspect_csv(id)
    try:
        db_connection.add_item_into_input(id, mechanism, nsp, systems)
        print(f"New input item registered: \"{id}\"")
    except Exception as e:
        print(e)
        print(f"Registration of input_id: \"{id}\" skipped")
    finally:
        return id

def take_results_files(execution_dir: str) -> list[str]:
    """Returns a list with the paths of the results files

    Args:
        results_dir (str): path of the directory where ODEApplication has used as output
        in the configuration yaml file

    Returns:
        list[str]: list with the paths of the results files
    """
    results = []
    for file in Path(execution_dir).iterdir():
        if file.is_dir():
            for result in file.iterdir():
                results.append(str(result))
            return results


if __name__ == "__main__":
    load_dotenv()

    # Runtime variables
    config_file: str = os.environ.get(CONFIG_YAML)
    workspace_dir: str = os.environ.get(EXECUTION_WORKSPACE)
    binary_app: str = os.environ.get(ODE_APPLICATION_BINARY)
    database_path: str = os.environ.get(DATABASE)

    # Database connection
    db = DataBase(database_path)

    execution_dir = prepare_execution_dir(workspace_dir)
    config_file_prepared = prepare_config_file(config_file, workspace_dir)
    input_id = register_input_if_needed(config_file_prepared, db)

    # Run integrator
    result = subprocess.run([binary_app, config_file_prepared], capture_output=True, text=True, check=True)
    print("ODEApplication output: " + result.stdout.strip("\n"))
    result_files = take_results_files(execution_dir)

    # TODO: Process files and save them into database

    db.close()
