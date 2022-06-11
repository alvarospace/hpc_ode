import numpy as np

class Model:
    def __init__(self, name, file_name) -> None:
        self.name = name
        self.file_name = file_name

    def read_data(self) -> None:
        self.data = np.loadtxt(self.file_name, delimiter=',', skiprows=1)



if __name__ == "__main__":
    ref_object = Model("ref", "ref_results.csv")
    ref_object.read_data()
    print(ref_object.data)