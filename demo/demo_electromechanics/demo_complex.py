""" 
A more complex example using the solver function (splitting solver) in the m3h3 package in 3D. 
"""
from m3h3 import *

from m3h3.pde.solver.electro_solver import BasicMonodomainSolver
from m3h3.setup_parameters import Physics, Parameters
from geometry import Geometry2D, Geometry

import matplotlib.pyplot as plt 

from cbcbeat import plot, Expression, Constant, BoxMesh, Point, refine, as_tensor, CompiledSubDomain, MeshFunction, Markerwise
from cbcbeat.utils import TimeStepper

import numpy as np

# Defining the mesh. 
Lx = 2. 
Ly = 7.
Lz = 3. 

dx = 0.1

N = lambda v: int(np.rint(v))
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(Lx, Ly, Lz), N(Lx/dx), N(Ly/dx), N(Lz/dx))
geo = Geometry(mesh)

C_m = 0.01
chi = 140

# Conductivities as defined by page 4339 of Niederer benchmark
sigma_il = 0.17  # mS / mm
sigma_it = 0.019 # mS / mm
sigma_el = 0.62  # mS / mm
sigma_et = 0.24  # mS / mm

# Compute monodomain approximation by taking harmonic mean in each
# direction of intracellular and extracellular part
def harmonic_mean(a, b):
    return a*b/(a + b)
sigma_l = harmonic_mean(sigma_il, sigma_el)
sigma_t = harmonic_mean(sigma_it, sigma_et)

# Scale conducitivites by 1/(C_m * chi)
s_l = sigma_l/(C_m*chi) # mm^2 / ms
s_t = sigma_t/(C_m*chi) # mm^2 / ms

# Define conductivity tensor
M = as_tensor(((s_l, 0, 0), (0, s_t, 0), (0, 0, s_t)))

time = Constant(0.0)

S1_marker = 1
L = 1.5
S1_subdomain = CompiledSubDomain("x[0] <= L + DOLFIN_EPS && x[1] <= L + DOLFIN_EPS && x[2] <= L + DOLFIN_EPS", L=L)
S1_markers = MeshFunction("size_t", mesh, mesh.topology().dim())
S1_subdomain.mark(S1_markers, S1_marker)

# Define stimulation (NB: region of interest carried by the mesh
# and assumptions in cbcbeat)
duration = 2. # ms
A = 50000. # mu A/cm^3
cm2mm = 10.
factor = 1.0/(chi*C_m) # NB: cbcbeat convention
amplitude = factor*A*(1./cm2mm)**3 # mV/ms
I_s = Expression("time >= start ? (time <= (duration + start) ? amplitude : 0.0) : 0.0",
                    time=time,
                    start=0.0,
                    duration=duration,
                    amplitude=amplitude,
                    degree=0)
# Store input parameters in cardiac model
stimulus = Markerwise((I_s,), (1,), S1_markers)

# Set up dt, starting time, and max time. 
dt = 0.1
t_0 = 0
t_max = 1
num_steps = int((t_max - t_0)/dt)
interval = (t_0, t_max)

# Define the conductivity (tensors). 
M_i = M
M_e = None

# Set up the parameteres for the model. 
params = Parameters("M3H3")
params.set_electro_parameters()
params.add("I_s" ,1)
params.add("M_i", M_i)
params.add("M_e", M_e) 
params.add("cell_model", "Tentusscher_panfilov_2006_M_cell")
params.add("time", Constant(t_0))
params.add("dt", dt)
params.add("I_a", 1)

# system = M3H3(geo, params, stimulus)

# # Run the simulation and print out the progress. 
# for (timestep, fields) in system.solve(interval, dt):
#     print("(t_0, t_1) = (%g, %g)", timestep)

#     (vs_, vs, vur) = fields


# # Plot the resulting solution fields. 
# plt.figure()
# plot(vs[0], title="Transmembrane potential (v) at end time")
# plt.savefig("TransmembranePot.png")
# plt.figure()
# plot(vs[-1], title="1st state variable (s_0) at end time")
# plt.savefig("s_0(T).png")

# print("Done!!")