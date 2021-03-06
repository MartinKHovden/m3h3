************************************************************************************
An introduction to doing cardiac simulations of cardiac electrophysiologoly in m3h3
************************************************************************************

This guide is a work in progress. The reader is assumed to have some knowledge
on how to use FEniCS and cbcbeat. For an introduction to FEniCS, see: 
https://fenicsproject.org/tutorial/. For an introudction to cbcbeat, see the 
documentation at: https://cbcbeat.readthedocs.io/en/latest/index.html. 

In this guide we will go through the most important steps in setting up an electrophysiologoly 
simulation in M3H3. This solver solves the monodomain or the bidomain equations presented in 
Sundnes (2006) for the transmembrance potentail along with state variables in the 
different cell models. 

Importing m3h3 
===============
To import all the functionality of m3h3, run 

.. code-block:: python

    from m3h3 import *

This will also import all the full functionality of FEniCS and cbcbeat.

Setting up the mesh
======================
The first step is to set up the computational domain of interest. This 
is done by creating a mesh using the functionality of FEniCS. All fenics 
mesh functions are included in the initial import of m3h3. 

Setting up a simple mesh can easily be done by using the functionality from 
FEniCS: 

.. code-block:: python 

    mesh = UnitSquareMesh(100, 100)

Alternatively, the mesh can be read from file. This can be done using the 
Mesh function.

.. code-block:: python 

    mesh = Mesh("data/mesh115_refined.xml.gz")

Other mesh-functions can be found in the mesh module of fenics: https://fenicsproject.org/docs/dolfin/2016.2.0/python/programmers-reference/cpp/mesh/index.html.

When the mesh it created, you need to set up a geometry object with 
this mesh

.. code-block:: python 

    geo = Geometry(mesh)

This will be necesarry when considering interactions in m3h3. 

Setting up the parameters and specifying the problems 
=======================================================
The Parameter class is where the user can specify what parameters they want to 
use in the simulations. It is a subclass of FEniCS' Parameter class and contains 
most of the information about the electrophysiologoly problem (and later also Fluid 
parameters, Solid parameters, and porous parameters). 

The Parameter class can be used like this 

.. code-block:: python 

    parameters = Parameters("M3H3")

When initiating the Parameters class it automatically sets the start-time and 
end-time to default values of 0 and 1. This can be changed as follows

.. code-block:: python 

    parameters["end_time"] = Constant(10.0)
    parameters["start_time"] = Constant(0.0)

to set the timespan of the simulations. start and end time is assumed to of type df.Constant. 
To add a nested set of electro parameters, 
you can run 

.. code-block:: python 

    parameters.set_electro_parameters()

This will set the electro parameters to the default values, which can be changed
in the same way as for the start- and end-time. Since no argument is given, 
the parameters are set to the default values. For more information on the function, 
see the API. 

We can now print out the parameters in the electro parameters to see what 
can be updated. 

.. code-block:: python  

    electro_params = parameters["Electro"]

Printing this to the terminal should give 

.. code-block:: python 

    ['stimulus', 'applied_current', 'initial_conditions', 'I_a', 'M_e', 'M_i', 'cell_model', 'dt', 'linear_variational_solver', 'pde_model', 'polynomial_degree', 'theta', 'use_average_u_constraint']

where the values can be updated as before. We will look closer at how to set up the
stimulus and applied current later in the guide. 

When setting the electro parameters, we can also change the solver parameters. Those 
have their own parameter set within the main parameter object that can be accessed:

.. code-block:: python 

    electrosolver_parameters = parameters["ElectroSolver"]


and again, printing out the keys gives

.. code-block:: python 

    ['BasicCardiacODESolver', 'BidomainSolver', 'CardiacODESolver', 'MonodomainSolver', 'apply_stimulus_current_to_pde', 'enable_adjoint', 'ode_solver_choice', 'pde_solver', 'theta']


Now that the electro parameters are set to default values, they can be changed 
as one would do in a dictionary. 

.. code-block:: python 

    electro_params["M_i"] = M_i
    electro_params["M_e"] = M_e
    electro_params["cell_model"] = "Beeler_reuter_1977"
    electro_params["dt"] = dt
    electro_params["stimulus"] = Expression("10*x[1]*t", t = Constant(0.0), degree = 1)

Note how the stimulus can be added to parameter set. We will look closer at how 
to set up the stimulus in a later section. 


For the electrical simulations a stimulus, applied current, and initial conditions
can be given. 

Now we can also change the parameters for the electro solver. 
This is done in a similar way as for the electro parameters 

.. code-block:: python 

    electrosolver_parameters["theta"] = 0.5                        # Second order splitting scheme
    electrosolver_parameters["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
    electrosolver_parameters["CardiacODESolver"]["scheme"] = "RL1" # 1st order Rush-Larsen for the ODEs
    electrosolver_parameters["MonodomainSolver"]["linear_solver_type"] = "iterative"
    electrosolver_parameters["MonodomainSolver"]["algorithm"] = "cg"
    electrosolver_parameters["MonodomainSolver"]["preconditioner"] = "sor"#"petsc_amg"
    electrosolver_parameters["apply_stimulus_current_to_pde"] = True

Here we see that we can choose between the monodomain and the bidomain equations. 

Stimulus 
++++++++++
The stimulus can be added as either a Constant, Expression, Markerwise function or a CompiledExpression. 
By using a Markerwise function or CompiledExpression, the position of the stimulus can be given. For more 
info on how to use subdomains and set up stimulus, see the FEniCS tutorial. In general,
all methods that work for setting up stimulus in cbcbeat also work in m3h3.   

Two examples of stimulus is shown below. The first is a simple stimulus using the 
Expression class. 

.. code-block:: python 

    stimulus = Expression("x[0]*t", t = start_time, degree = 1)

This is a simple stimulus that moves along the x-axis with time. Note that 
the time is set equal to the start_time object. This is so that the stimulus 
time is syncronized with the internal time of m3h3.  

A more complex example uses the CompiledSubdomain functionality in combination 
with the Markerwise class to set up two separate stimulus in the domain.
The first step is to mark the two areas of the domain where the stimulus should 
be applied. 

.. code-block:: python

    stimulus_domain = MeshFunction("size_t", mesh, mesh.topology().dim())
    stimulus_domain.set_all(0)

    stimulus_1 = CompiledSubDomain("pow(x[0],2) + pow(x[1],2) <= 0.5 + tol", tol = 1e-15 )
    stimulus_1.mark(stimulus_domain, 1)

    stimulus_2 = CompiledSubDomain("pow(x[0]-1.0,2) + pow(x[1]-1, 2) <= 0.1 + tol", tol = 1e-15)
    stimulus_2.mark(stimulus_domain, 2)

When the two subdomains are set up, the stimulus for each domain can be set the following way 

.. code-block:: python 

    I_s_1 = Expression("t >= start ? (t <= (duration + start) ? amplitude : 0.0) : 0.0",
                t=start_time,
                start=0.0,
                duration=1,
                amplitude=10,
                degree=0)

    I_s_2 = Expression("t >= start ? (t <= (duration + start) ? amplitude : 0.0) : 0.0",
                t=start_time,
                start=0.0,
                duration=0.5,
                amplitude=5,
                degree=0)

Note that the string in expression can be any expression allowed in c++.
Also note that the time variable in the expression is set to the start_time object. 
This is important so that the time is syncronized with the internal time of m3h3.

The stimuluses can now be applied to the subdomains via the Markerwise class 

.. code-block:: python 

    stimulus = Markerwise((I_s_1, I_s_2), (1,2), stimulus_domain)

Alternativelly, it is possible to use the CompiledExpression function from fenics. 
To see an example on how this is done, see the demo folder. 

Setting up m3h3 
==================
Now that all the parameters are set, we can create an instance of the 
m3h3 class. 

.. code-block:: python 

    system = m3h3(geo, params)

Running the simulation 
=======================
The m3h3 object can now be used to run the simulations. There are two different 
ways of doing this: Using the step method, or using the solve method in m3h3.  

Running simulations with the step function 
++++++++++++++++++++++++++++++++++++++++++++
To run the simulations using the step function, we have to know the number of 
steps to run. In the parameter object, the start and end time is stored, as 
well as the step length. The number of steps can then be calculated

.. code-block:: python 

    num_steps = int((end_time - start_time)/dt)

To run the simulations, we can set up a for loop that runs the step function 
for each iteration

.. code-block:: python 

    for _ in range(num_steps):
        print("Time interval: ", (float(system.time), float(system.time) + dt) )
        system.step()


This will also print out the time interval it solves for for each iteration. 
Each call to the step function updates the solution fields of system. Those can 
be extracted using the get_solution_field() function

.. code-block:: python 

    vs_, vs, vur = system.get_solution_field()["Electro"]

where we are only interested in the solution fields for the electro problem. 


Running simulations with the solve function 
+++++++++++++++++++++++++++++++++++++++++++++
Alternatively, it is possible to use the solve function for doing the same 
simulation. The solve function calls the step function multiple times. It returns 
a generator that can be iterated over to obtain the solution fields and 
the time intervals. 

.. code-block:: python 
    
    for (t0, t1), solution_field in system.solve():
        print((t0, t1))
        vs_, vs, vur = solution_field

Again, the solution fields can be extracted using the get_solution_field() function
as we did for the step function. 

Post-processing 
================
The last part is to visualize the results. There are different ways of doing this, 
and it depends on the problem dimension. 

Plotting in 2D
+++++++++++++++++
The easiest way to plot when looking at a 2D problem is to use the plot function
from fenics. The plot function depends on matplotlib. If you dont have 
matplotlib installed on your system, it can easily be obtained by using pip 

.. code-block:: python 

    pip install matplotlib 

To plot the results, you can run the plot function with the desired field to 
plot 

.. code-block:: python 

    plot(vs[0], title="Plot of transmembrane potential")

This will plot the transmenbrane potential over the domain. 

Plotting in 3D
+++++++++++++++++
In the previous example where the mesh was taken from a file, the domain 
is in 3 dimensions. The easiest way to visualize the results in 3D is to use 
external plotting software. Two of the possibilities is to 
use ParaView or vedo. 

To download ParaView, follow the instructions on: https://www.paraview.org/download/
When the plotting software is installed on the system, we need a file to 
visualize. FEniCS have a function called File() that can convert the solution fields 
into various formats. For visualization in ParaView and vedo, .pvd(tvk) files are 
a possible format that can be used. To write the output to file, use the 
File() function from FEniCS

.. code-block:: python 

    File("filename.pvd") << vs.split()[0]

filename.pvd can now be found in the present folder and opened using ParaView. 

If you want a lighter package to do similar plotting using the terminal, you can use the vedo 
python package. This can be installed using pip 

.. code-block:: python

    pip install -U vedo 

Then you can open a python script and run 

.. code-block:: python 

    from vedo import show 
    show("./filename.pvd")

This assumes that you are in the folder where the output data is stored. 