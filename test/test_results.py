import numpy as np
import os
import sys
from tabulate import tabulate

class Model:
    # Constructor
    def __init__(self, file_name):
        self.file_name = file_name
        self.name = os.path.splitext(file_name)[0].split('_')[-1]
        
        # Read data
        self.data = np.loadtxt(self.file_name, delimiter=',', skiprows=1)
        self.data = self.data[:,0:53]


    # Calculation of the Frobenius norm
    def diff_norm(self, model):
        return np.linalg.norm(self.data - model.data, 'fro')



def print_table(rows, header):
    pass

if __name__ == "__main__":

    # Check for input reference file
    ref_filename = "results_reference.csv"
    if len(sys.argv) > 1:
        ref_filename = sys.argv[1]

    # Read folder files and take ".csv"
    files = os.listdir()
    results = []
    for file in files:
        if file.endswith(".csv"):
            results.append(file)
    
    # Select reference file and remove from target list
    
    results.remove(ref_filename)

    # Output Table
    header = ["Files", "Frobenius Norm"]
    rows = []

    # Process data
    ref_model = Model(ref_filename)
    for file in results:
        target = Model(file)
        res = ref_model.diff_norm(target)
        rows.append([target.name, res])

    print("\nReference file:   ", ref_filename)
    print("\nFrobenius Norm calculation to validation:\n")

    # Print table in console
    print(tabulate(rows, headers=header))

