import cantera as ct
import numpy as np

# Add name of mechanism file 
gas = ct.Solution('../mechanisms/ESTiMatE-Mech_04042022.xml')

# Create a database based on Temperature and Amount of fuel
# Array of Temperature 
T = np.linspace(1000,2100,100)
# Array of Fuel Concentration
phi = np.linspace(0.5,3.5,100)

# Output file name 
f = open('res_est.csv','w')

# Write Header 
for nn in range(gas.n_species):
  f.write("CO{},".format( nn ))
f.write("ENTHA,")
f.write("TEMPE\n")

# Write Actual Data 
for i in T:
  for j in phi:
    gas.set_equivalence_ratio(j,'N-C12H26:1','N2:3.76,O2:1')
    gas.TP = i, ct.one_atm
    for nn in range(gas.n_species):
# Write Concentrations 
      f.write("{},".format( gas.Y[nn] ))
# Write Enthalpy
    f.write("{},".format( gas.enthalpy_mass ))
# Write Temperature
    f.write("{}".format( gas.T ))
    f.write("\n")
f.close()
