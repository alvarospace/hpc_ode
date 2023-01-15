import os
import subprocess
import shutil
from pathlib import Path

import yaml
from dotenv import load_dotenv

from common import get_repo_directory, inspect_csv
from database import DataBase
from constants import (
    INPUT_TABLE,
    EXECUTION_DIR,
    CONFIG_YAML,
    ODE_APPLICATION_BINARY
)

DATA = "ref_data/res_gri3.0.csv"

def prepare_execution_folder(directory: str) -> None:
    """Remove from filesystem what is in "directory" and prepare
    the directory for the execution

    Args:
        directory (str): path to the directory where the program will use for
        the temporary and result files in the execution
    """
    dir_path = Path(directory)
    if dir_path.exists():
        shutil.rmtree(str(dir_path))
    dir_path.mkdir(parents=True)
    out_path = dir_path.joinpath("out")
    out_path.mkdir()


if __name__ == "__main__":
    load_dotenv()
    # input_path = Path(get_repo_directory()).joinpath(DATA)
    # nsp, systems = inspect_csv(str(input_path))
    # print(nsp, systems)
    # db = DataBase("./ODEIntegratorDB.db")
    # data = db.query(f"SELECT * FROM {INPUT_TABLE}")
    # print(data)
    # db.close()

    # result = subprocess.run(["echo", "Hola desde python"], capture_output=True, text=True, check=True)
    # print(result.stdout.strip("\n"))

    prepare_execution_folder(EXECUTION_DIR)
    config_file: str = os.environ.get(CONFIG_YAML)

    # Change outFolder name of the config_file
    with open(config_file, "r") as f:
        yaml_dict: dict = yaml.safe_load(f)
    yaml_dict["outFileService"]["outFolder"] = EXECUTION_DIR
    # Write changes to the 
    yaml.dump

    print(config_file)
    binary: str = os.environ.get(ODE_APPLICATION_BINARY)
    print(binary)

    # result = subprocess.run([binary, config_file], capture_output=True, text=True, check=True)
    # print(result.stdout.strip("\n"))


