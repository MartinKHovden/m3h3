from dolfin import (Constant)

from m3h3 import Physics
from m3h3.pde.solver import Solver


class FluidSolver(Solver):

    def __init__(self, form, time, interval, dt, parameters, **kwargs):
        super().__init__(form, time, interval, dt, parameters, **kwargs)
