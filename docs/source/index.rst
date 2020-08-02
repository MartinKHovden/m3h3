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
the heart. This documentation and the framework is a work in progress. 

M3H3 builds on the functionality of software specific for solving the 
various parts that goes into the simulation of the heart. This is the 
electrical activity, fluid flow, continuum mechanics and ... 

CBCBeat is used for simulations of the electrical activity. See 
https://cbcbeat.readthedocs.io/en/latest/index.html for more info. The 
framework is based on the agorithms described in Sundnes et al (2006). 

The software includes functionality for doing full simulations of the 
electrical activity, fluid flow, and solid properties of the heart. 

The installation page shows how to obtain m3h3. At the bottom of this page
a simple example of how to use the framework is shown. To get more details
on how the framework can be used, see the guide and demos. 

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
