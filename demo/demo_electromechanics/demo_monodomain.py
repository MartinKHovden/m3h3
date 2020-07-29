""" 
A simple example on how to use the m3h3 module for solving for the electrical 
activity of the heart. 

This example shows how to:
- Set up a simple geometry 
- Update the parameters of the electro problem
- Update the solver parameters for the electro solver 
- Add a simple stimulus
- Run the electro simulation 
- Plot the results 

"""

from m3h3 import *
import matplotlib.pyplot as plt 

# Define the computational domain.
mesh = UnitSquareMesh(100, 100)

# Set up the geometry given the computational domain.
geo = Geometry2D(mesh)

# Set up dt, t_0, and t_max: 
dt = 0.1
t_0 = 0.0
t_max = 1.0
num_steps = int((t_max - t_0)/dt)
interval = (t_0, t_max)

# Define the conductivity (tensors):
M_i = 1.0
M_e = 1.0

# Set up the parameteres for the heart-model: 
params = Parameters("M3H3")

# Set the end and start time of the simulation: 
params["end_time"] = t_max
params["start_time"] = t_0 

# Set the parameters for the electro problem: 
params.set_electro_parameters()
electro_params = params["Electro"]
electro_params["M_i"] = M_i
electro_params["M_e"] = M_e
electro_params["cell_model"] = "Beeler_reuter_1977"#"Tentusscher_panfilov_2006_M_cell"
electro_params["dt"] = dt
electro_params["stimulus"] = Expression("10*x[1]*t", t = Constant(0.0), degree = 1)

# Set the electro solver parameters: 
electrosolver_params = params["ElectroSolver"]
electrosolver_params["theta"] = 0.5                        # Second order splitting scheme
electrosolver_params["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
electrosolver_params["CardiacODESolver"]["scheme"] = "RL1" # 1st order Rush-Larsen for the ODEs
electrosolver_params["MonodomainSolver"]["linear_solver_type"] = "iterative"
electrosolver_params["MonodomainSolver"]["algorithm"] = "cg"
electrosolver_params["MonodomainSolver"]["preconditioner"] = "petsc_amg"

# Initialize the system with parameters and geometry:
system = M3H3(geo, params)

# Run the simulation by using the step function:
for i in range(num_steps):
    print("Time interval: ", (float(system.time), float(system.time) + dt) )

    system.step()

# Extract the solution:
vs_, vs = system.get_solution_fields()[str(Physics.ELECTRO)]

# Plot the resulting solution fields:
plt.figure()
plot(vs[0], title="Transmembrane potential (v) at end time")
plt.savefig("TransmembranePot.png")
plt.figure()
plot(vs[-1], title="1st state variable (s_0) at end time")
plt.savefig("s_0(T).png")

print("Done!!")