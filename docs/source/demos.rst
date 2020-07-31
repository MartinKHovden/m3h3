Demos
==========
This guide provides an basic introduction on how to run a simple simluation 
of the electrical activity in the heart using the M3H3 framework. The normal 
workflow when doing simulation with M3H3 is 

1. Import M3H3 as 

    .. code-block:: python
        :linenos:

        from m3h3 import *

2. Set up the mesh and the geometry
    .. code-block:: python 
        :linenos:

        mesh = UnitSquareMesh(100, 100)
        geo = Geometry2D(mesh)

3. 
