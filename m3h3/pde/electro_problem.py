# -*- coding: utf-8 -*-
"""This module keeps track of the cell models and the cardiac models for the 
electro problem. It also stores the solution fields for the electro problem. 
"""

import cbcbeat
from cbcbeat.cardiacmodels import CardiacModel

from m3h3.pde import Problem
import m3h3.ode

import os
import importlib

class ElectroProblem(Problem):
    """This class sets up the electro problem with cell model and 
    cardiac model as well as storing the solution fields. 

    """

    def __init__(self, geometry, time, parameters, **kwargs):
        super().__init__(geometry, time, parameters, **kwargs)

        # Set up the stimulus for the electro problem:
        self.stimulus = None
        self._set_up_stimulus(**kwargs)

        # Set up the applied current for the electro problem: 
        self.applied_current = None
        self._set_up_current(**kwargs)

        # Set up the initial conditions for the electro problem: 
        self.cell_model_initial_conditions = None
        self._set_up_initial_conditions(**kwargs)

        # Not sure how to handle the cell model params: 
        self.cell_model_parameters = None 

        # Sets up the cell model and cardiac model. 
        self.cell_model = self.get_cell_model()
        self.cardiac_model = self.get_cardiac_model()



    def get_cell_model(self):
        """Returns the cell model specified in the parameters.
        
        Iterates over all files in the ode-folder (which is where the cell 
        models are stored) and checks if the cell model specified in 
        the parameteres are implemented. If not, it raises an error and returns
        a list of the implemented cell models. All implemented cell models 
        should also be included in __all__ in __init__.py in the ode-folder.

        Assumes that all cell-models are implemented in the ode-folder. The 
        naming convention is that the .py-file is named the same as the 
        cell-class in the file, but with a lower-case first letter. The 
        class name is the same but with an upper-case first letter.  

        """

        # Get the cell model name from the user parameters: 
        model = self.parameters['cell_model']

        cell_model = None

        # Iterates over all files in the ode-folder to see if model is implemented:
        for filename in os.listdir(str(os.path.dirname(__file__)) + "/../ode"):
            # Checks if filename is a Python-file: 
            if filename.endswith(".py"):
                filename_split = filename.split(".")
                # Checks if model name is same as filename:
                if filename_split[0].lower() == model.lower():
                    # Get the cell model from the corresponding file: 
                    cell_model = getattr(importlib.import_module("m3h3.ode." + 
                                        filename_split[0]), 
                                        filename_split[0][:1].capitalize() +
                                        filename_split[0][1:])
                            
        # If the cell model is found in the ode-folder, return it:
        if cell_model != None:
            return cell_model(params = self.cell_model_parameters, 
                    init_conditions = self.cell_model_initial_conditions)
        # If cell model not in ode-folder, raise NotImplementedError: 
        else: 
            raise NotImplementedError("""Cell model not implemented, 
                        try one of the following:""", m3h3.ode.__all__)

    def get_cardiac_model(self):
        """Returns the cardiac model for the electro problem given the cell 
        model and the user-parameters.  

        """  
        print(self.parameters["M_i"], self.parameters["M_e"])
        print(self.stimulus)
        print(self.time)
        return CardiacModel(domain = self.geometry.mesh,
                                time = self.time, 
                                M_i = self.parameters["M_i"], 
                                M_e = self.parameters["M_e"], 
                                cell_models = self.cell_model, 
                                stimulus = self.stimulus,
                                applied_current = self.applied_current)

    # def update_solution_fields(self, solution, prev_current):
    def update_solution_fields(self, vs_, vs, vur):
        """Function for updating the solution field. Used in step 
        function of m3h3. 
        """
        self.vs_ = vs_ 
        self.vs = vs
        self.vur = vur


    def _get_solution_fields(self):
        """ Return the solution field for the electro problem.

        Function for returning the solution fields for 
        
        """
        return (self.vs_, self.vs, self.vur)


    def _set_up_stimulus(self, **kwargs):
        """Add the given stimulus to the electro problem. Stimulus is 
        suposed to be of type CompiledExpression, Expression or Markerwise.

        *Note*
            When using CompiledExpression, time should be encoded as t in the 
            expression. When using Expression or Markerwise, the time can be 
            encoded as t or time in the user_parameters of the expression or 
            markerwise. The time variable in the expression is then set to 
            the internal time of the m3h3 object. 

        """ 
        self.stimulus = self.parameters["stimulus"]

        # if self.stimulus != None:
        #     if isinstance(self.stimulus, cbcbeat.Markerwise):
        #         for stim in self.stimulus.values():
        #             if "t" in stim.user_parameters:
        #                 stim.t = self.time
        #             elif "time" in stim.user_parameters:
        #                 stim.time = self.time

        #     elif isinstance(self.stimulus, cbcbeat.Constant):
        #         pass

        #     elif isinstance(self.stimulus, cbcbeat.Expression):
        #         if "t" in self.stimulus.user_parameters:
        #             self.stimulus.t = self.time
        #         elif "time" in self.stimulus.user_parameters:
        #             self.stimulus.time = self.time

        #     elif isinstance(self.stimulus, cbcbeat.CompiledExpression):
        #         self.stimulus.t = self.time._cpp_object

        #     else:
        #         msg = """Stimulus should be an Expression, CompiledExpression 
        #         or Markerwise, not %r""" %type(self.stimulus)
        #         raise TypeError(msg)


    def _set_up_current(self, **kwargs):
        """ Add the given applied current to the electro problem. 

        *Note*
            Time variable should be encoded as time or t. 
            The time variable is then set to the internal time of the 
            m3h3 solver. 

        """
        self.applied_current = self.parameters["applied_current"]

        # if "applied_current" in kwargs.keys():
        if self.applied_current != None:
            if "t" in self.applied_current.user_parameters:
                self.applied_current.t = self.time
            elif "time" in self.applied_current.user_parameters:
                self.applied_current.time = self.time

    def _set_up_initial_conditions(self, **kwargs):
        """ Set the initial conditions for the problem. 
        """

        self.cell_model_initial_conditions = self.parameters["initial_conditions"]


