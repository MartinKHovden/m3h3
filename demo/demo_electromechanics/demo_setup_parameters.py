""" This demo shows how to use the parameter class in m3h3. 

"""

from m3h3 import *
# from dolfin import Parameters

# Initiate the parameter object with the label "M3H3"
params = Parameters("M3H3")

# Check which parameters are in params:
print("keys in params: ", params.keys())
# Results in: ['end_time', 'log_level', 'start_time']

# To see the value for each key:
print("end_time: ", params["end_time"])
# Results in: 1.0

# Those can be changed in the following way:
params["end_time"] = 10.0

# To see the updated value:
print("end_time after update: ", params["end_time"])
# Results in: 10.0

# To add more parameters, use
params.add("test", 10.0)

# and print out the keys to see the set update
print("keys after adding test:", params.keys())

# To add a parameter set for the electro parameters and the electro solver parameters:
params.set_electro_parameters()

# The electro parameter set and the electro solver set can be accesed this way:
electro_parameters = params["Electro"]
electro_solver_parameters = params["ElectroSolver"]

print("1", params.keys())
params.add("tst", 1)
print(params.keys())

# The default values can then be 

print(electro_parameters.keys())
print(electro_solver_parameters.keys())



