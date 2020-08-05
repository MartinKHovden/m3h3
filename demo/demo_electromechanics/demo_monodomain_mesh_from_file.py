""" 
An example on how to use m3h3 to do simulations of the electrical 
activity of the heart for a more complex stimulus applied to the domain.
The mesh is read from file.  

This example shows how to:
- Read the geometry from file. 
- Update the parameters of the electro problem.
- Update the solver parameters for the electro solver.
- Add a more complex stimulususing Meshfunctions and Compiled subdomains. 
- Run the electro simulation.
- Plot the results.

"""

from cbcbeat import *
from m3h3 import *

import matplotlib.pyplot as plt 

mesh = Mesh("data/mesh115_refined.xml.gz")
mesh.coordinates()[:] /= 1000.0 # Scale mesh from micrometer to millimeter
mesh.coordinates()[:] /= 10.0   # Scale mesh from millimeter to centimeter
mesh.coordinates()[:] /= 4.0    # Scale mesh as indicated by Johan/Molly

geo = Geometry(mesh)

# Set up the stimulus subdomains: 
stimulation_cells = MeshFunction("size_t", mesh, "data/stimulation_cells.xml.gz")

time = Constant(0.0)

V = FunctionSpace(mesh, "DG", 0)
from stimulation import cpp_stimulus
amp = 30 # stimulus amplitude
pulse = CompiledExpression(compile_cpp_code(cpp_stimulus).Stimulus(),
                            element=V.ufl_element(), t=time._cpp_object,
                            amplitude=amp, duration=10.0,
                            cell_data=stimulation_cells)

# Set up dt, t_0, and t_max: 
dt = 0.1
t_0 = 0.0
t_max = 1.0
num_steps = int((t_max - t_0)/dt)
interval = (t_0, t_max)

# Define the conductivity (tensors):
M_i = 2.0
M_e = 1.0

# Set up the parameteres for the heart-model: 
params = Parameters("M3H3")

print(params.keys())

params["end_time"] = t_max
params["start_time"] = t_0 

params.set_electro_parameters()

electro_params = params["Electro"]
electro_params["dt"] = dt
electro_params["M_i"] = M_i
electro_params["M_e"] = M_e
electro_params["cell_model"]  = "Beeler_reuter_1977"#"Tentusscher_panfilov_2006_M_cell"
electro_params["stimulus"]= pulse
electro_params["applied_current"] = None

# Set up the parameters for the splitting solver: 
electrosolver_parameters = params["ElectroSolver"]
electrosolver_parameters["theta"] = 0.5                        # Second order splitting scheme
electrosolver_parameters["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
electrosolver_parameters["CardiacODESolver"]["scheme"] = "RL1" # 1st order Rush-Larsen for the ODEs
electrosolver_parameters["MonodomainSolver"]["linear_solver_type"] = "iterative"
electrosolver_parameters["MonodomainSolver"]["algorithm"] = "cg"
electrosolver_parameters["MonodomainSolver"]["preconditioner"] = "sor"#"petsc_amg"
electrosolver_parameters["apply_stimulus_current_to_pde"] = True

# Initialize the system with parameters and geometry.
system = M3H3(geo, params, time = time)

# Run the simulation by using the step function:
for i in range(num_steps):
    print("Time interval: ", (float(system.time), float(system.time) + dt) )
    system.step()

# Or run the simulations by using the solve function: 
# for (t0, t1), solution_field in system.solve():
#     print((t0, t1))

# Extract the solution:
vs_, vs = system.get_solution_fields()[str(Physics.ELECTRO)]

File("test.pvd") << vs.split()[0]

# The results can be visualized using vedo or ParaView. 
# If you are running on wsl, you migh have to install vedo in a 
# windows terminal and then show the file from there. 

print("Done!!")