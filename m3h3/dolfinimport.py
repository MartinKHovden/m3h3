# FIXME: This is here for readthedocs Mock purposes. Better fix would
# be, duh, better.
from dolfin import (Parameters, Mesh, Constant, Expression, assemble,
                    LUSolver, KrylovSolver, PETScKrylovSolver,
                    dx, Measure, parameters, VectorFunctionSpace, Function,
                    DirichletBC, TrialFunction, TestFunction, solve, inner,
                    warning, error, begin, end, set_log_level)

from ufl.log import info_blue, info_red, info_green, debug
from ufl import zero
from dolfin.cpp.function import GenericFunction

try:
    from dolfin import UserExpression, CompiledExpression
except:
    pass

try:
    from dolfin.cpp.log import log, LogLevel
except:
    pass
    
from dolfin import *
import dolfin

try:
    from dolfin_adjoint import *
    import dolfin_adjoint

except:
    # FIXME: Should we raise some sort of warning?
    dolfin_adjoint = None
    pass
