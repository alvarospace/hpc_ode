from pathlib import Path

COMMON_RELATIVE_LOCATION = "/app/common.py"

def get_repo_directory() -> str:
    """Returns repository absolute path

    Returns:
        str: Absolute path of the hpc_cvode repository in the current host
    """
    return __file__.replace(COMMON_RELATIVE_LOCATION, "")

def inspect_csv(filename: str) -> tuple[int, int]:
    """Inspect csv file and returns the number of species
    and the number of systems that it contains, it does
    not return the name of the mechanism

    Args:
        filename (str): Path to the csv file

    Returns:
        tuple[int, int]: First item is the number of species and second
        item is the number of systems
    """
    with open(filename, "r") as csvfile:
        data = csvfile.readlines()
        header: list[str] = data[0].split(",")
    
    systems = len(data) - 1
    header_cleaned = [item.strip('"') for item in header]
    nsp = 0
    for item in header_cleaned:
        if item.startswith("CO"):
            nsp += 1
    return nsp, systems

def get_mechanism_list() -> list[str]:
    """Returns a list of the available mechanisms

    Returns:
        list[str]: List of strings with the .yaml mechanism files available
    """
    resources_path = Path(get_repo_directory()).joinpath("resources/mechanisms")
    available_mechanisms: list[str]= []
    for file in resources_path.iterdir():
        if file.suffix == ".yaml":
            available_mechanisms.append(file.name)

    return available_mechanisms
