""" 
A simple example for solving the Monodomain equations using 
a splitting solver from the m3h3 package.  
"""

from m3h3 import *

from m3h3.pde.solver.electro_solver import BasicMonodomainSolver
from m3h3.setup_parameters import Physics, Parameters
from geometry import Geometry2D

import matplotlib.pyplot as plt 

from cbcbeat import plot, Expression, Constant
from cbcbeat.utils import TimeStepper

# Define the computational domain
mesh = UnitSquareMesh(100, 100)
time = Constant(0.0)

# Set up the geometry given the computational domain 
geo = Geometry2D(mesh)

# Set up time step, starting time, and ending time 
dt = 0.1
t_0 = 0
t_max = 1
num_steps = int((t_max - t_0)/dt)
interval = (t_0, t_max)

# Define the conductivity (tensors)
M_i = 2.0
M_e = 1.0

# Set up the parameteres for the model 
params = Parameters("M3H3")
params.set_electro_parameters()
params.add("I_s" ,1)
params.add("M_i", M_i)
params.add("M_e", M_e)
params.add("cell_model", "Tentusscher_panfilov_2006_M_cell")
params.add("time", Constant(t_0))
params.add("dt", dt)
params.add("I_a", 1)
params.add("stimulus", Expression("10*t*x[0]", t=time, degree=1))

# Initialize the system 
system = M3H3(geo, params)

# t = t_0

for (timestep, fields) in system.solve(interval, dt):
    print("(t_0, t_1) = (%g, %g)", timestep)

    # Extract the components of the field
    (vs_, vs, vur) = fields

# Visualize some results
plt.figure()
plot(vs[0], title="Transmembrane potential (v) at end time")
plt.savefig("TransmembranePot.png")
plt.figure()
plot(vs[-1], title="1st state variable (s_0) at end time")
plt.savefig("s_0(T).png")

print("Done!!")