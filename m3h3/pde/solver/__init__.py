from m3h3.pde.solver.solver import Solver
# from m3h3.pde.solver.electro_solver import BasicBidomainSolver
# from cbcbeat import SplittingSolver
from m3h3.pde.solver.electro_solver import ElectroSolver
from m3h3.pde.solver.solid_solver import SolidSolver
from m3h3.pde.solver.fluid_solver import FluidSolver
from m3h3.pde.solver.porous_solver import PorousSolver

__all__ = ['Solver', 'ElectroSolver', 'SolidSolver', 'FluidSolver',
            'PorousSolver']