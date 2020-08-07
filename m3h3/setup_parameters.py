# -*- coding: utf-8 -*-
"""This module handles the parameters for cardiac simulations.

It contains one main class for storing the general parameters for each problem.
In addition, it contains a class for the parameters of the electro problem. This
is because df.Parameters only take floats, ints and strings. Since we 
might want to add complex stimuluses, an alternative is needed. By using the 
ElectroParameter class, this is now possible. 

For more information on how to use the parameter class, see the user guide in 
the documentaion at m3h3-docs.readthedocs.com. 

"""

from enum import Enum

import dolfin as df
from dolfin import (LogLevel, LUSolver, PETScKrylovSolver,
                    NonlinearVariationalSolver)

import cbcbeat

import m3h3


class Physics(Enum):
    """This Enum contains physics descriptors for cardiac simulations.
    """
    ELECTRO = "Electro"
    SOLID = "Solid"
    FLUID = "Fluid"
    POROUS = "Porous"

    def __str__(self):
        return self.value

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

def set_dolfin_compiler_parameters():
    """Sets dolfin parameters to speed up the compiler.
    """
    flags = ["-O3", "-ffast-math", "-march=native"]
    df.parameters["form_compiler"]["quadrature_degree"] = 4
    df.parameters["form_compiler"]["representation"] = "uflacs"
    df.parameters["form_compiler"]["cpp_optimize"] = True
    df.parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)

class ElectroParameters(df.Parameters):
    """Class for handling the electro parameters.

    This class handles the electro parameters. It inherits most of 
    the functionality from df.Parameters, but is modified to be able 
    to store more complex stimulus, applied current, and conductivities. 
    This is needed because df.Parameters only takes floats, ints or strings. 
    Instead of storing stimulus, applied current, and conductivities in the 
    actual df.Parameter storage, it stores them in separate variables. The 
    class is updated so that those items can be accessed in the same way 
    as items in a normal df.Parameter object. 

    *Arguments*
        label (:py:class:`str`)
            The label of the parameter set.  

    *Notes*
        When creating an instance of ElectroParameters, the stimulus, applied 
        current, and initial conditions are already in the set with default
        values of None. Those can updated as usual.  

    *Usage*
        .. code-block::python

            params = ElectroParameters("ElectroParams")      
    """
    def __init__(self, label, **kwargs):
        super().__init__(label, **kwargs)

        # Add storage for the problem specifications:
        self.stimulus = None 
        self.applied_current = None
        self.initial_conditions = None
        self.M_i = None
        self.M_e = None
    
    def add(self, *args):
        """ Function for adding new parameters. Should only be used 
        for adding new elements to the set. For updating already defined 
        specifications, use:

        .. code-block:: python 

            problem_specifications["specification"] = Something...

        """
        if args[0] == "stimulus": 
            print("Stimulus is already in the parameter set") 
            raise ValueError("Stimulus is already in the parameter set")
        elif args[0] == "applied_current":
            print("Applied current is already in the parameter set")
            raise ValueError("Applied current is already in the parameter set")
        elif args[0] == "initial_conditions":
            print("Initial conditions is already in the parameter set")
            raise ValueError("Initial condition is already in the parameter set")
        elif args[0] == "M_i":
            print("M_i is already in the parameter set")
            raise ValueError("M_i is already in the parameter set")
        elif args[0] == "M_e":
            print("M_e is already in the parameter set")
            raise ValueError("M_e is already in the parameter set")
        else :
            super().add(*args)

    def __getitem__(self, key):
        if key == "stimulus":
            return self.stimulus 
        elif key == "applied_current":
            return self.applied_current
        elif key == "initial_conditions":
            return self.initial_conditions
        elif key == "M_i":
            return self.M_i 
        elif key == "M_e":
            return self.M_e
        else:
            return super().__getitem__(key)
    
    def __setitem__(self, key, value):
        if "stimulus" == key:
            self.stimulus = value
        elif key == "applied_current":
            self.applied_current = value
        elif key == "initial_conditions":
            self.initial_conditions = value
        elif key == "M_i":
            self.M_i = value
        elif key == "M_e":
            self.M_e = value
        else:
            super().__setitem__(key, value)

    def get(self, key):
        """Function for getting the value for a key.

        *Returns*
            Returns the value corrensponding to the key

        """
        if key == "stimulus":
            return self.stimulus 
        elif key == "applied_current":
            return self.applied_current
        elif key == "initial_conditions":
            return self.initial_conditions
        elif key == "M_i":
            return self.M_i 
        elif key == "M_e":
            return self.M_e
        else :
            return super().get(key)    

    def keys(self):
        """Returns a dictionary with the keys in the parameter set. 
        
        *Returns*
            return_list (:py:class:`dict`)
                List with the keys in the parameter set. 
        """
        return_list = ["stimulus", "applied_current", "initial_conditions", "M_i", "M_e"]
        return_list += super().keys()
        
        return return_list

    def has_key(self, key):
        """Checks if key is in parameter set. 

        Function for checking if key is in the parameter set. Returns True if it 
        is in the parameter set and false otherwise. 

        *Returns*
            (:py:class:`boolean`)
                Boolean that represents if key is in parameter set or not. 

        """
        if key == "stimulus":
            return True
        elif key == "applied_current":
            return True
        elif key == "initial_conditions":
            return True
        elif key == "M_i":
            return True
        elif key == "M_e":
            return True
        else :
            return super().has_key(key)

    def clear(self):
        """Clears the parameter set.

        Function for removing every parameter from the parameter set. 

        *Note*
            The special parameters (stimulus, applied current, M_i, M_e) are 
            set to None. 

        """
        self.stimulus = None 
        self.applied_current = None 
        self.initial_conditions = None
        self.M_i = None
        self.M_e = None
        super().clear()

class Parameters(df.Parameters):
    """Class for handling the parameters for the cardiac simulations. 
    
    This class handles the parameters for the cardiac simulations. It 
    inherits most of the functionality from  `dolfin`'s Parameters class. 
    The main difference is that this class have an alternative way 
    of handling the electro parameters, since we should be able to add more complex
    stimuluses (df.Parmeters does only allows floats, ints and strings).
    Therefore, the ElectroParameter class is used, and are stored in a separate
    variable self.electro_parameters instead of in the actual df.Parameter dict. 

    The parameter set can be nested. Each problem have their own nested
    Parameter set that contains the parameters. 

    *Notes*
        start_time and end_time is assumed to be of type df.Constant. 

    *Usage*

        The Parameters object can be initialized with a given name 

        .. code-block:: python 

            params = Parameters("M3H3")

        To see the keys in the parameter set, use the key() function:

        .. code-block:: python 

            keys = params.keys()


        Parameters can be added as follows:

        .. code-block:: python 

            params.add("param1", 1.0)

        They can be changed as follows:

        .. code-block:: python 

            params["param1"] = 2.0

        They can be extracted as follows:

        .. code-block:: python 

            p1 = params["param1"]

        To set the default parameters for the electro problem, do

        .. code-block:: python 

            params.set_electro_parameters()

        This will add a new parameter set with electro parmeters to params. 
        To change the parameters in the electro parameter set, do 

        .. code-block:: python 

            electro_params = params["Electro"]
        
        The electro parameters can then be viewed, added, changed, and extracted as before. 

    """
    def __init__(self, label, **kwargs):
        """Initializes a new instance of the Parameter class.

        Calls the super-class initializer, and automatically sets 
        the dolfin compiler parameters as well as the default parameters.
        See self.set_default_parameters for more info. Also creates a 
        new variable to store the electro parameter set. This is necessary since
        the electro paramters needs its own class. See :py:class:`ElectroParamters`
        for more info. 

        """
        super().__init__(label, **kwargs) 
        set_dolfin_compiler_parameters()
        self.set_default_parameters()

        self.electro_parameters = None
        self.start_time = None
        self.end_time = None

    def set_default_parameters(self):
        """Sets default simulation parameters.

        Function for setting the default parameters. Default parameters are 
        "start_time", "end_time" and "log_level". This is automatically 
        called when initializing an object of m3h3. 

        *Note*
            This is automatically called when initializing a new instance, so 
            does not need to be used. Also note that start_time and end_time
            is assumed to be of type df.Constant. 
        """
        self.add("log_level", df.get_log_level())
        m3h3.log(self["log_level"], "Log level is set to {}".format(
            LogLevel(self["log_level"])
        ))
        self.add("start_time", df.Constant(0.0))
        self.add("end_time", df.Constant(1.0))

    def __getitem__(self, key):
        if key == Physics.ELECTRO.value:
            return self.electro_parameters
        elif key == "start_time":
            return self.start_time 
        elif key == "end_time":
            return self.end_time
        else:
            return super().__getitem__(key)

    def __setitem__(self, key, value):
        if key == Physics.ELECTRO.value:
            self.electro_parameters = value
        elif key == "start_time":
            self.start_time = value 
        elif key == "end_time":
            self.end_time = value
        else:
            super().__setitem__(key, value)

    def keys(self):
        """Returns the keys in the parameter set.

        *Returns*
            :py:class:`dict`
                Dictionary containing all the keys in the parameter set. 
        """
        keys = super().keys()
        if self.electro_parameters != None:
            print("Here")
            keys = keys + [Physics.ELECTRO.value]
        if self.start_time != None:
            keys = keys + ["start_time"]
        if self.end_time != None:
            keys = keys + ["end_time"]
        return keys

    def has_parameter_set(self, parameter_set):
        """ Checks if parameter set is in Parameters object. 

        *Arguments*
            :py:class:`string`
                Name of parameter set to check if in object. 
        """ 
        if (parameter_set == Physics.ELECTRO.value and
                            self.electro_parameters != None):
            return True
        elif (parameter_set == Physics.ELECTRO.value and
                            self.electro_parameters == None):
            return False 
        else:
            return super().has_parameter_set(parameter_set)

    def set_electro_parameters(self, parameters=None):
        """Sets parameters for electrophysiology problems and solver. If
        argument is None, default parameters are applied.

        """
        if not self.has_parameter_set(Physics.ELECTRO.value):
            self._set_electro_default_parameters()
            self._set_electro_solver_default_parameters()
        if self.has_parameter_set(Physics.ELECTRO.value) and parameters:
            self[Physics.ELECTRO.value].update(parameters)

    def set_solid_parameters(self, parameters=None):
        """Sets parameters for solid mechanics problems and solver. If
        argument is None, default parameters are applied.
        """
        if not self.has_parameter_set(Physics.SOLID.value):
            self._set_solid_default_parameters()
        if self.has_parameter_set(Physics.SOLID.value) and parameters:
            self[Physics.SOLID.value].update(parameters)

    def set_fluid_parameters(self, parameters=None):
        """Sets parameters for fluid mechanics problems and solver. If
        argument is None, default parameters are applied.
        """
        if not self.has_parameter_set(Physics.FLUID.value):
            self._set_fluid_default_parameters()
        if self.has_parameter_set(Physics.FLUID.value) and parameters:
            self[Physics.FLUID.value].update(parameters)

    def set_porous_parameters(self, parameters=None):
        """Sets parameters for porous mechanics problems and solver. If
        argument is None, default parameters are applied.
        """
        if not self.has_parameter_set(Physics.POROUS.value):
            self._set_porous_default_parameters()
        if self.has_parameter_set(Physics.POROUS.value) and parameters:
            self[Physics.POROUS.value].update(parameters)

    def _set_electro_default_parameters(self):
        """Sets the default parameters for the electro problem.
        """
        self.electro_parameters = ElectroParameters(Physics.ELECTRO.value)

        # Set default parameters
        self.electro_parameters.add("dt", 0.001)
        self.electro_parameters.add("theta", 0.5)
        self.electro_parameters.add("polynomial_degree", 1)
        self.electro_parameters.add("use_average_u_constraint", False)
        self.electro_parameters.add("I_a", 0.0)
        self.electro_parameters.add("cell_model", "Tentusscher_panfilov_2006_M_cell")
        self.electro_parameters.add("pde_model", "bidomain")

        self.electro_parameters["M_i"] = 1.0
        self.electro_parameters["M_e"] = 2.0

        self.electro_parameters.add(df.LinearVariationalSolver.default_parameters())

    def _set_electro_solver_default_parameters(self):
        """ Sets the default splitting solver parameters to be used in the
        simulation. Default parameters are given by the splittingsolver class 
        in cbcbeat. 
        """ 
        electro_solver = cbcbeat.splittingsolver.SplittingSolver.default_parameters()
        electro_solver.rename("ElectroSolver")
        self.add(electro_solver)

    def _set_solid_default_parameters(self):
        solid = df.Parameters(Physics.SOLID.value)
        solid.add("dt", 1e-3)
        solid.add("dummy_parameter", False)

        # Add boundary condtion parameters
        solid.add(df.Parameters("BoundaryConditions"))
        solid["BoundaryConditions"].add("base_bc", "fixed")
        solid["BoundaryConditions"].add("lv_pressure", 10.0)
        solid["BoundaryConditions"].add("rv_pressure", 0.0)
        solid["BoundaryConditions"].add("pericardium_spring", 0.0)
        solid["BoundaryConditions"].add("base_spring", 0.0)

        # Add default parameters from both LU and Krylov solvers
        solid.add(NonlinearVariationalSolver.default_parameters())
        solid.add(LUSolver.default_parameters())
        solid.add(PETScKrylovSolver.default_parameters())

        # Add solver parameters
        solid.add(df.Parameters("Solver"))
        solid["Solver"].add("dummy_parameter", False)

        self.add(solid)

    def _set_fluid_default_parameters(self):
        fluid = df.Parameters(Physics.FLUID.value)
        fluid.add("dt", 1e-3)
        fluid.add("dummy_parameter", False)

        # Add default parameters from both LU and Krylov solvers
        fluid.add(LUSolver.default_parameters())
        fluid.add(PETScKrylovSolver.default_parameters())

        # Add solver parameters
        fluid.add(df.Parameters("Solver"))
        fluid["Solver"].add("dummy_parameter", False)

        self.add(fluid)

    def _set_porous_default_parameters(self):
        porous = df.Parameters(Physics.POROUS.value)
        porous.add("dt", 1e-3)
        porous.add("dummy_parameter", False)

        # Add default parameters from both LU and Krylov solvers
        porous.add(LUSolver.default_parameters())
        porous.add(PETScKrylovSolver.default_parameters())

        # Add solver parameters
        porous.add(df.Parameters("Solver"))
        porous["Solver"].add("dummy_parameter", False)

        self.add(porous)
