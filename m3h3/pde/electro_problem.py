# -*- coding: utf-8 -*-
"""This module implements the variational form for electrophysiology problems
"""

from dolfin import (grad, inner, Constant, FiniteElement, Function,
                    FunctionAssigner, FunctionSpace, Measure, MixedElement,
                    TestFunctions, TrialFunctions, UserExpression, Expression)
import cbcbeat
from cbcbeat.cardiacmodels import CardiacModel

from m3h3.pde import Problem
# from m3h3.ode import Tentusscher_panfilov_2006_M_cell
import m3h3.ode

import os
import importlib


class Stimulus(UserExpression):

    def __init__(self, markers, stimulus_marker, **kwargs):
        super().__init__(degree=kwargs['degree'])
        self.markers = markers
        self.stimulus_marker = stimulus_marker
        self.amplitude = kwargs['amplitude']
        self.period = kwargs['period']
        self.duration = kwargs['duration']
        self.t = kwargs['t']

    def eval_cell(self, values, x, cell):
        periodic_t = float(self.t) % self.period
        # if self.markers[cell.index] == self.stimulus_marker\
        #                                     and periodic_t < self.duration:
        #     values[0] = self.amplitude
        if self.markers["STIMULUS"] == self.stimulus_marker\
                                            and periodic_t < self.duration:
            values[0] = self.amplitude
        else:
            values[0] = 0
            
    def value_shape(self):
        return ()


class ElectroProblem(Problem):
    """This class implements the variational form for electrophysiology
    problems.

    It initiates the cardiac cell model and the cardiac model with the 
    given stimulus and applied current. 
    """

    def __init__(self, geometry, time, parameters, **kwargs):
        super().__init__(geometry, time, parameters, **kwargs)

        self.stimulus = None
        self._set_up_stimulus(**kwargs)

        self.applied_current = None
        self._set_up_current(**kwargs)

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

        """
        model = self.parameters['cell_model']

        cell_model = None

        # FIXME: Check if everything works properly when adding new files. 
        # Also try to make it more robust if user provides different 
        # capitalization of words. 


        # Iterates over all files in the ode-folder where the cell models are 
        # stored:
        for filename in os.listdir(str(os.path.dirname(__file__)) + "/../ode"):
            if filename.endswith(".py"):
                filename_split = filename.split(".")
                if filename_split[0].lower() == model.lower():
                    cell_model = getattr(importlib.import_module("m3h3.ode." + 
                                        filename_split[0]), 
                                        filename_split[0][:1].capitalize() +
                                        filename_split[0][1:])
                            
        if cell_model != None:
            return cell_model()
        else: 
            raise NotImplementedError("""Cell model not implemented, 
                        try one of the following:""", m3h3.ode.__all__)


    def get_cardiac_model(self):
        """Returns the cardiac model for the electro problem given the cell 
        model and parameters.  
        """ 
        print(self.parameters["M_i"])
        return CardiacModel(domain = self.geometry.mesh,
                                time = self.time, 
                                M_i = self.parameters["M_i"], 
                                M_e = self.parameters["M_e"], 
                                cell_models = self.cell_model, 
                                stimulus = self.stimulus,
                                applied_current = self.applied_current)


    def update_solution_fields(self, solution, prev_current):
        """ Function for updating the solution field. Used in step 
        function of m3h3. 
        """
        self.solution = solution 
        self.prev_current = prev_current


    def _get_solution_fields(self):
        return (self.prev_current, self.solution)


    def _set_up_stimulus(self, **kwargs):
        """ Add the given stimulus to the electro problem. Stimulus is 
        suposed to be of type Expression or Markerwise. 
        """ 
        if "stimulus" in kwargs.keys():
            if isinstance(kwargs["stimulus"], cbcbeat.Markerwise):
                self.stimulus = kwargs["stimulus"]
                for stim in self.stimulus.values():
                    if "t" in stim.user_parameters:
                        stim.t = self.time
                    elif "time" in self.applied_current.user_parameters:
                        self.applied_current.time = self.time

            elif isinstance(kwargs["stimulus"], cbcbeat.Expression):
                self.stimulus = kwargs["stimulus"]
                if "t" in self.stimulus.user_parameters:
                    self.stimulus.t = self.time
                elif "time" in self.applied_current.user_parameters:
                    self.applied_current.time = self.time

            else:
                msg = """Stimulus should be an Expression 
                or Markerwise, not %r""" %type(kwargs["stimulus"])
                raise TypeError(msg)


    def _set_up_current(self, **kwargs):
        """ Add the given applied current to the electro problem. Applied 
        current is suposed to be of type ufl.Expr. 
        """
        if "applied_current" in kwargs.keys():
            self.applied_current = kwargs["applied_current"]
            if "t" in self.applied_current.user_parameters:
                self.applied_current.t = self.time
            elif "time" in self.applied_current.user_parameters:
                self.applied_current.time = self.time


