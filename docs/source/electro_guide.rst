*********************************************************************************
An introduction to doing cardiac simulations of the electrical activity in M3H3
*********************************************************************************

This guide is a work in progress. 

In this guide we will go through the most important steps in setting up a cardiac
simulation in M3H3. Currently, only the electrical activity can be simulated. 
This solver solves the monodomain or the bidomain equations presented in 
Sundnes (2006) to find the transmembrane potential. 

The goal is for the framework to be able to combine simulations of the 
electrical activity with the simulations of the solid mechanics, fluid 
mechanics, as well as the porous properties of the heart. For each step in 
setting up the simulation, we will describe it from the most basic to more 
advanced. For full working examples, see the demos folder. 

Importing M3H3 
===============
To import all the functionality of M3H3, run 

.. code-block:: python

    from m3h3 import *

This will also import the full functionality of FEniCS and CBCBeat, so those 
do not have to be imported separately. 

Setting up the mesh
======================
The first step is to set up the computational domain of interest. This 
is done by creating a mesh using the functionality of fenics. All fenics 
mesh functions are included in the initial import of m3h3. 

Setting up a simple mesh can easily be done by using the functionality from 
feincs: 

.. code-block:: python 

    mesh = UnitSquareMesh(100, 100)

Alternatively, the mesh can be read from file. This can be done using the 
Mesh function.

.. code-block:: python 

    mesh = Mesh("data/mesh115_refined.xml.gz")


Setting up the parameters and specifying the problems 
=======================================================
The Parameter class is where the user can specify what parameters they want to 
use in the simulations. It is a subclass of fenics' Parameter class and contains 
most of the information about the electro problem (and later also Fluid 
parameters, Solid parameters, and porous parameters). The Parameter class 
contains nested Parameter classes that contains the parameters for each 
of the problems. The Parameter class can be used like this 

.. code-block:: python 

    parameters = Parameters("M3H3")

When initiating the Parameters class it automatically sets the start-time and 
end-time to default values of 0 and 1. This can be changed as follows

.. code-block:: python 

    parameters["end_time"] = 10.0
    parameters["start_time"] = 0.0

to decide the timespan of the simulations. To add a nested set of electro parameters, 
you can run 

.. code-block:: python 

    parameters.set_electro_parameters()

This will set the electro parameters to the default values, which can be changed
in the same way as for the start- and end-time. Since no argument is given, 
the parameters are set to the default values. For more information on the function, 
see the API. 

Setting up the electro simulations
++++++++++++++++++++++++++++++++++++++++

Now that the parametres contains a nested set of electro parameters, they can 
be changed as before. 

.. code-block:: python 

    electro_parameters = params["Electro"]

There are multiple ways to set the parameters for the electro simulations. The 
easiest is to first set them equal to the default electro parameters and then 
updating them from there.

Now that the electro parameters are set to default values, they can be changed 
as one would do in a dictionary. 

.. code-block:: python 

    electro_params = params["Electro"]
    electro_params["M_i"] = M_i
    electro_params["M_e"] = M_e
    electro_params["cell_model"] = "Beeler_reuter_1977"
    electro_params["dt"] = dt
    electro_params["stimulus"] = Expression("10*x[1]*t", t = Constant(0.0), degree = 1)

Note how the stimulus can be added to parameter set. We will look closer at how 
to set up the stimulus in a later section. 

Alternatively, they can be set manually as a function argument: 

.. code-block:: python 

    dfvdn


For the electrical simulations a stimulus, applied current, and initial conditions
can be given. 

Now we can also change the parameters for the electro solver. 
This is done in a similar way by 



Stimulus 
++++++++++
The stimulus can be added as either a Constant, Expression, or a Markerwise function. 
By using a Markerwise function, the position of the stimulus can be given. For more 
info on how to use subdomains and set up stimuluses, see the FEniCS tutorial.  

Two examples of stimulus is shown below. The first is a simple stimulus using the 
Expression class. 

.. code-block:: python 

    stimulus = Expression("x[0]*t", t = Constant(0.0), degree = 1)

This is a simple stimulus that moves along the x-axis with time.  

A more complex example uses the CompiledSubdomain functionality in combination 
with the Markerwise class to set up two separate stimuluses in the domain.
The first step is to mark the two areas of the domain where the stimuluses should 
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
                t=Constant(0.0),
                start=0.0,
                duration=1,
                amplitude=10,
                degree=0)

    I_s_2 = Expression("t >= start ? (t <= (duration + start) ? amplitude : 0.0) : 0.0",
                t=Constant(0.0),
                start=0.0,
                duration=0.5,
                amplitude=5,
                degree=0)

Note that the string in expression can be any expression allowed in c++. The stimuluses can now be connected to the subdomains via the Markerwise class 

.. code-block:: python 

    stimulus = Markerwise((I_s_1, I_s_2), (1,2), stimulus_domain)



Setting up the fluid simulations 
+++++++++++++++++++++++++++++++++++

Setting up the porous simulations 
+++++++++++++++++++++++++++++++++++++

Setting up the interactions
++++++++++++++++++++++++++++++

Running the simulation 
=======================

Post-processing 
================
The last part is to 