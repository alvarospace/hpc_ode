from pathlib import Path

from common import get_repo_directory
from database import DataBase
from constants import INPUT_TABLE

DATA = "ref_data/res_gri3.0.csv"

def inspect_csv(filename: str) -> tuple[str, str]:
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

if __name__ == "__main__":
    input_path = Path(get_repo_directory()).joinpath(DATA)
    nsp, systems = inspect_csv(str(input_path))
    print(nsp, systems)


