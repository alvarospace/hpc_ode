import os
import subprocess
import shutil
import datetime
from pathlib import Path

import yaml
from dotenv import load_dotenv

from common import inspect_csv
from database import DataBase
from runtime_variables import (
    EXECUTION_WORKSPACE,
    CONFIG_YAML,
    ODE_APPLICATION_BINARY,
    DATABASE
)

PERFORMANCE_FILE = "performance.yaml"
CONFIG_FILE = "config.yaml"
LOG_FILE = "out.log"
RESULT_FILE = ".csv"

def read_yaml(filename: str) -> dict:
    with open(filename, "r") as f:
        yaml_dict: dict = yaml.safe_load(f)
    return yaml_dict

def to_binary(filename: str) -> bytes:
    with open(filename, "rb") as f:
        bin_file = f.read()
    return bin_file

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
    
    print("Preparing configuration file...")
    # Change outFolder name of the config_file
    yaml_dict = read_yaml(config_filename)
    out_destination_path = Path(destination_path).joinpath("out")
    yaml_dict["outFileService"]["outFolder"] = str(out_destination_path)

    # Write changes to the destination path
    config_tmp_path = Path(destination_path).joinpath(CONFIG_FILE)
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
    config_yaml = read_yaml(config_filename)
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

def take_result_files(execution_dir: str) -> dict[str, str]:
    """Returns a list with the paths of the results files

    Args:
        results_dir (str): path of the directory where ODEApplication has used as output
        in the configuration yaml file

    Returns:
        dict[str, str]: dict with the paths of the results files
    """
    results = {}
    for file in Path(execution_dir).iterdir():
        if file.is_dir():
            for result in file.iterdir():
                name = result.name
                if name == CONFIG_FILE:
                    results[CONFIG_FILE] = result
                elif name == PERFORMANCE_FILE:
                    results[PERFORMANCE_FILE] = result
                elif name == LOG_FILE:
                    results[LOG_FILE] = result
                elif name.endswith(RESULT_FILE):
                    results[RESULT_FILE] = result
                else:
                    print("No match taking the result files")
            return results

def main():
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
    result_files = take_result_files(execution_dir)

    # Dictionaries to insert into database
    config_yaml = read_yaml(result_files[CONFIG_FILE])
    performance_yaml = read_yaml(result_files[PERFORMANCE_FILE])
    
    execution_dict = {
        "input_id": input_id,
        "reader": config_yaml["reader"]["type"],
        "writer": config_yaml["writer"]["type"],
        "integrator": config_yaml["integrator"]["type"],
        "logger": config_yaml["logger"]["type"],
        "output": to_binary(result_files[RESULT_FILE]),
        "read_time": performance_yaml["readTime"],
        "integration_time": performance_yaml["integrationTime"],
        "write_time": performance_yaml["writeTime"],
        "total_time": performance_yaml["totalTime"],
    }

    integrator_dict = {
        "reltol": config_yaml["integrator"]["reltol"],
        "abstol": config_yaml["integrator"]["abstol"],
        "pressure": config_yaml["integrator"]["pressure"],
        "dt": config_yaml["integrator"]["dt"]
    }

    log_binary = result_files.get(LOG_FILE)
    if log_binary is not None:
        log_binary = to_binary(log_binary)
    log_dict = {
        "log_level": config_yaml["logger"]["logLevel"],
        "log_file": log_binary
    }

    omp_dict: dict | None = None
    if str(config_yaml["integrator"]["type"]).endswith("OMP"):
        omp_dict = {
            "cpus": config_yaml["integrator"]["omp"]["cpus"],
            "schedule": config_yaml["integrator"]["omp"]["schedule"]["type"],
            "chunk": config_yaml["integrator"]["omp"]["schedule"]["chunk"]
        }
    
    db.insert(execution_dict, integrator_dict, log_dict, omp_dict)
    db.close()

if __name__ == "__main__":
    main()
