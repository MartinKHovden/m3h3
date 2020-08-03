.. M3H3 documentation master file, created by
   sphinx-quickstart on Thu Jul  9 10:01:57 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#############################################
M3H3: A framework for cardiac simulations. 
#############################################

Welcome to M3H3's documentation. M3H3 is an Python-based open source framework for doing cardiac simulations. 
It build on the finite element solver FEniCS and combines functionality from 
various frameworks for simulating different processes in the heart. 
The goal is for M3H3 to be able to combine simulations of electrical 
activivty, solid mechanics, fluid mechanics, and the porous properties of 
the heart. This documentation and m3h3 is a work in progress. 

The documentation is split into different parts. In the design specification, 
the aim of the framework is described as reference for further development of the library. 
In the installation guide different methods of obtaining the software is described. 
In the user guide we go through the main parts in setting up various simulations 
in m3h3. In demos, full examples of running programs are shown. Finally, the 
API reference gives an detailed overview of the classes and functions in m3h3. 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   designspecifications
   installation
   user_guide
   demos
   API


Quickstart
------------------------------

An general example of how the M3H3 library can be used for doing cardiac 
simulation:

.. code-block:: python
   :linenos:

   from m3h3 import *

   mesh = UnitSquareMesh(100, 100)

   geo = Geometry2D(mesh)

   system = m3h3(geo) 

   for (time, solution_fields) in system.solve():
      # do something with the solution. 
   
   plot(solution_field[0])

This simple example creates a square domain of 100 points in each dimension. 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
