from constants import INPUT_DATA

COMMON_RELATIVE_LOCATION = "app/common.py"

def get_absolute_path_input_data() -> list:
    input = INPUT_DATA
    repo_directory = __file__.replace(COMMON_RELATIVE_LOCATION, "")
    for entry in input:
        entry["input_id"] = repo_directory + entry["input_id"]
    return input