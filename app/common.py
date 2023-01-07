from pathlib import Path

from constants import INPUT_DATA

COMMON_RELATIVE_LOCATION = "/app/common.py"

def get_repo_directory() -> str:
    return __file__.replace(COMMON_RELATIVE_LOCATION, "")

def get_absolute_path_input_data() -> list:
    input = INPUT_DATA
    repo_directory = __file__.replace(COMMON_RELATIVE_LOCATION, "")
    for entry in input:
        entry["input_id"] = repo_directory + entry["input_id"]
    return input

def get_mechanism_list() -> list[str]:
    resources_path = Path(get_repo_directory()).joinpath("resources/mechanisms")
    available_mechanisms: list[str]= []
    for file in resources_path.iterdir():
        if file.suffix == ".yaml":
            available_mechanisms.append(file.name)

    return available_mechanisms
