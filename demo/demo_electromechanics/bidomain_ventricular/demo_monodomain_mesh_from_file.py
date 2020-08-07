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

The set up for this problem is taken from the demos in cbcbeat, written by: Marie E. Rognes. 

"""

from cbcbeat import *
from m3h3 import *

import matplotlib.pyplot as plt 

# Read the mesh from file
mesh = Mesh("data/mesh115_refined.xml.gz")
mesh.coordinates()[:] /= 1000.0 # Scale mesh from micrometer to millimeter
mesh.coordinates()[:] /= 10.0   # Scale mesh from millimeter to centimeter
mesh.coordinates()[:] /= 4.0    # Scale mesh as indicated by Johan/Molly

geo = Geometry(mesh)

def set_up_conductivities(healthy = True):
    Vv = VectorFunctionSpace(mesh, "DG", 0)
    fiber = Function(Vv)
    File("data/fibers.xml.gz") >> fiber
    sheet = Function(Vv)
    File("data/sheet.xml.gz") >> sheet
    cross_sheet = Function(Vv)
    File("data/cross_sheet.xml.gz") >> cross_sheet

    # Extract stored conductivity data.
    V = FunctionSpace(mesh, "CG", 1)
    if (healthy == True):
        info_blue("Using healthy conductivities")
        g_el_field = Function(V, "data/healthy_g_el_field.xml.gz", name="g_el")
        g_et_field = Function(V, "data/healthy_g_et_field.xml.gz", name="g_et")
        g_en_field = Function(V, "data/healthy_g_en_field.xml.gz", name="g_en")
        g_il_field = Function(V, "data/healthy_g_il_field.xml.gz", name="g_il")
        g_it_field = Function(V, "data/healthy_g_it_field.xml.gz", name="g_it")
        g_in_field = Function(V, "data/healthy_g_in_field.xml.gz", name="g_in")
    else:
        info_blue("Using ischemic conductivities")
        g_el_field = Function(V, "data/g_el_field.xml.gz", name="g_el")
        g_et_field = Function(V, "data/g_et_field.xml.gz", name="g_et")
        g_en_field = Function(V, "data/g_en_field.xml.gz", name="g_en")
        g_il_field = Function(V, "data/g_il_field.xml.gz", name="g_il")
        g_it_field = Function(V, "data/g_it_field.xml.gz", name="g_it")
        g_in_field = Function(V, "data/g_in_field.xml.gz", name="g_in")

    # Construct conductivity tensors from directions and conductivity
    # values relative to that coordinate system
    A = as_matrix([[fiber[0], sheet[0], cross_sheet[0]],
                    [fiber[1], sheet[1], cross_sheet[1]],
                    [fiber[2], sheet[2], cross_sheet[2]]])
    from ufl import diag
    M_e_star = diag(as_vector([g_el_field, g_et_field, g_en_field]))
    M_i_star = diag(as_vector([g_il_field, g_it_field, g_in_field]))
    M_e = A*M_e_star*A.T
    M_i = A*M_i_star*A.T

    return M_i, M_e

# Set up dt, t_0, and t_max. 
dt = 0.1
start_time = Constant(0.0)
end_time = Constant(1.0)
num_steps = int((float(end_time) - float(start_time))/dt)
M_i, M_e = set_up_conductivities(healthy=True)

# Set up the stimulus subdomains from file: 
stimulation_cells = MeshFunction("size_t", mesh, "data/stimulation_cells.xml.gz")

V = FunctionSpace(mesh, "DG", 0)
from stimulation import cpp_stimulus
pulse = CompiledExpression(compile_cpp_code(cpp_stimulus).Stimulus(),
                            element=V.ufl_element(), t=start_time._cpp_object,
                            amplitude=30.0, duration=10.0,
                            cell_data=stimulation_cells)

# Set up the parameteres for the heart-model: 
params = Parameters("M3H3")

params["end_time"] = end_time
params["start_time"] = start_time 

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
electrosolver_parameters["pde_solver"] = "bidomain"  #"monodomain"        # Use Monodomain model for the PDEs
electrosolver_parameters["CardiacODESolver"]["scheme"] = "RL1" # 1st order Rush-Larsen for the ODEs
electrosolver_parameters["MonodomainSolver"]["linear_solver_type"] = "iterative"
electrosolver_parameters["MonodomainSolver"]["algorithm"] = "cg"
electrosolver_parameters["MonodomainSolver"]["preconditioner"] = "sor"#"petsc_amg"
electrosolver_parameters["apply_stimulus_current_to_pde"] = True

# Initialize the system with parameters and geometry.
system = M3H3(geo, params)

# Run the simulation by using the step function:
for i in range(num_steps):
    print("Time interval: (%.2f, %.2f)"% (float(system.time), float(system.time) + dt) )
    system.step()

# Extract the solution:
vs_, vs, vur = system.get_solution_fields()[str(Physics.ELECTRO)]

File("test.pvd") << vs.split()[0]

# The results can be visualized using vedo or ParaView. 
# If you are running on wsl, you migh have to install vedo in a 
# windows terminal and then show the file from there. 

print("Done!!")