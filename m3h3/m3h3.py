import numpy as np

from dolfin import (Constant, Parameters)

from geometry import HeartGeometry, MultiGeometry

from m3h3.setup_parameters import Parameters, Physics
from m3h3.pde import ElectroProblem, SolidProblem, PorousProblem, FluidProblem
from m3h3.pde.solver import *

from cbcbeat.cardiacmodels import CardiacModel
from cbcbeat import Expression, plot

""" 
TODO:
Fix how to cell models are searched through and given in params. 
Fix readthedocs, make the code look better. 
Check if stimulus works properly. Ok way to add stimulus.
Update initial conditions. 
Explain how to add more cell models (lower case in filename, upper case for class name.)
Write unit tests. 

M3H3 is a framework used for modelling and simulating the heart. It inherits 
much of the functionality from other libaries and combines it into a full 
framework that encompasses the different branches of physics relevant to 
caridac modelling. The goal is for the framework to be able to simulate 
the cardiac mechanics, cardiac electrophysiology, hemodynamic modelling, etc.

For now, only methods for simluating the electrical activity is implemented. 
For modelling the electrical activity in the heart, the splittingsolver from 
cbcbeat is used. Most of the functionality is inherited from the cbcbeat package.
It solves the coupled heart equations, presented in Sundnes et. al.:
https://www.researchgate.net/publication/265487224_Computing_the_Electrical_Activity_in_the_Human_Heart#:~:text=The%20contraction%20of%20the%20heart,called%20the%20electrocardiogram%20(ECG).

"""
class M3H3(object):
    """ A class for representing the full cardiac system. 

    To add stimulus, add it as a keyword argument with the keyword "stimulus". 
    If stimulus is dependent on time, use keyword "t" or "time" to connect it 
    to the internal timer in the m3h3-object. The stimulus should either be of 
    type Expression or Markewise. Use Markerwise if the position of the 
    stimulus should be used as well as if multiple stimuluses should be applied. 

    *Arguments*
        geometry :py:class:`geometry.Geometry`
            A geometry object that contains the mesh of the simulation domain. 
        parameters (:py:class:`dolfin.Parameters`)
            A Parameters object that contains parameters for setting up the 
            cardiac simulation. 
    

    The idea is to set up a M3H3 object that can be used to run simulations. 
    The basic usage is:
        - Set up the geometry with a mesh using the fenics.geometry package.  
        - Set the parameters for the cardiac model. 
        - Set the solver parameters. 
        - Create the M3H3 object given the geometry and parameters
        - Use M3H3 object's step/solve functions for running the simulations. 
        - Do post-processing of the output. 
    
    """
    def __init__(self, geometry, parameters, *args, **kwargs):
        self.parameters = parameters
        self.physics = [Physics(p) for p in parameters.keys()
                                                    if Physics.has_value(p)]
        self.interactions = kwargs.get('interactions', [])
        self.t_step_lengths = self._get_physics_dt()

        if len(self.interactions) > 0:
            self._check_physics_interactions()

        if 'time' in kwargs.keys():
            self.time = kwargs['time']
            kwargs.pop('time', None)
        else:
            self.time = Constant(self.parameters['start_time'])

        self._setup_geometries(geometry, self.physics)
        self._setup_problems(**kwargs)
        self._setup_solvers(**kwargs)



    def step(self):
        """ Does one step for the solvers. Finds the interval using the internal 
        time variable. If each problems solver uses different size of time step, 
        it runs the number of steps to make up for the difference. 

        Usage:
        system = m3h3(...)
        for i in range(num_steps):
            pre-process the data. 
            system.step()
            post-process the data.
        """

        # Setup the number of steps for each solver if it is the first iteration. 
        if self.time.values()[0] == self.parameters["start_time"]:
            self.num_step, self.max_dt = self._get_num_steps()

        time = float(self.time)

        dt = self._get_physics_dt()
        
        if Physics.ELECTRO in self.physics:
            for _ in range(self.num_step[str(Physics.ELECTRO)]):
                # Interval to solve for: 
                interval = (time, time+dt[str(Physics.ELECTRO)])

                # Does one step and extracts the solution fields. 
                self.electro_solver.step(interval)
                solution_fields = self.electro_solver.solution_fields()

                # Updates m3h3's solution fields. 
                self.electro_problem.update_solution_fields(solution_fields[0],
                                                            solution_fields[1])  
        
        if Physics.SOLID in self.physics:
            pass
    
        if Physics.FLUID in self.physics:
            pass 

        if Physics.POROUS in self.physics:
            pass 
        
        self.time.assign(time + self.max_dt)


    def solve(self):
        """ Solves the problem by running multiple steps for the full interval 
            ("start_time", "end_time"), where "start_time" and "end_time" are 
            given as parameters to the m3h3 object. 

            Returns a generator that can be iterated over:
            system = m3h3(...)
            for (timestep, fields) in system.solve():
                do something with the fields. 

            # FIXME: Update so that it uses the internal step function instead. 
        """
        max_dt = self._get_num_steps()
        end_time = self.parameters["end_time"]
        t0 = self.time
        t1 = t0 + max_dt

        while self.time + max_dt <= end_time:
            self.step()
            yield (t0, t1), self.get_solution_fields()
            t0 = t1 
            t1 += max_dt


    def get_solution_fields(self):
        """ Returns the solution fields for the different problems. Return a 
        dictionary with keys given by strings representing which problem they 
        correspond to. 
        """ 

        solution_fields = {}
        if Physics.ELECTRO in self.physics:
            solution_fields[str(Physics.ELECTRO)] =\
                                    self.electro_problem._get_solution_fields()
        if Physics.SOLID in self.physics: 
            pass
        if Physics.FLUID in self.physics:
            pass
        if Physics.POROUS in self.physics:
            pass
        return solution_fields


    def _get_num_steps(self):
        dt_physics = self._get_physics_dt()
        min_dt = min(dt_physics.values())
        max_dt = max(dt_physics.values())
        num_steps = {}
        if Physics.ELECTRO in self.physics:
            dt = dt_physics[Physics.ELECTRO]
            if self._check_dt_is_multiple(dt, min_dt):
                num_steps[Physics.ELECTRO] = int(max_dt/dt)
        if Physics.SOLID in self.physics:
            dt = dt_physics[Physics.SOLID]
            if self._check_dt_is_multiple(dt, min_dt):
                num_steps[Physics.SOLID] = int(max_dt/dt)
        if Physics.FLUID in self.physics:
            dt = dt_physics[Physics.FLUID]
            if self._check_dt_is_multiple(dt, min_dt):
                num_steps[Physics.FLUID] = int(max_dt/dt)
        if Physics.POROUS in self.physics:
            dt = dt_physics[Physics.POROUS]
            if self._check_dt_is_multiple(dt, min_dt):
                num_steps[Physics.POROUS] = int(max_dt/dt)
        return num_steps, max_dt


    def _check_dt_is_multiple(self, dt, min_dt):
        """ Function for checking if two time steps for different 
        problems are multiples of each other. 
        """ 
        if not (dt/min_dt).is_integer():
            msg = "Time step sizes have to be multiples of each other."\
                    "{} is not a multiple of {}".format(dt, min_dt)
            raise ValueError(msg)
        else:
            return True


    def _get_physics_dt(self):
        """ Returns a dictionary containing the time steps for each 
        problems solver. 
        """ 
        dt = {}
        if Physics.ELECTRO in self.physics:
            dt[Physics.ELECTRO] = self.parameters[str(Physics.ELECTRO)]['dt']
        if Physics.SOLID in self.physics:
            dt[Physics.SOLID] = self.parameters[str(Physics.SOLID)]['dt']
        if Physics.FLUID in self.physics:
            dt[Physics.FLUID] = self.parameters[str(Physics.FLUID)]['dt']
        if Physics.POROUS in self.physics:
            dt[Physics.POROUS] = self.parameters[str(Physics.POROUS)]['dt']
        return dt


    def _setup_problems(self, **kwargs):
        """ Set up the problems with the given parameters and geometry. 
        """
        if Physics.ELECTRO in self.physics:
            self.electro_problem = ElectroProblem(
                                        self.geometries[Physics.ELECTRO],
                                        self.time,
                                        self.parameters[str(Physics.ELECTRO)],
                                        **kwargs, problem_specifications = 
                                        self.parameters["problem_specifications"])
            
        if Physics.SOLID in self.physics:
            self.solid_problem = SolidProblem(
                                        self.geometries[Physics.SOLID], 
                                        self.time,
                                        self.parameters[str(Physics.SOLID)],
                                        **kwargs)

        if Physics.FLUID in self.physics:
            self.fluid_problem = FluidProblem(
                                        self.geometries[Physics.FLUID],
                                        self.time,
                                        self.parameters[str(Physics.FLUID)],
                                        **kwargs)

        if Physics.POROUS in self.physics:
            self.porous_problem = PorousProblem(
                                        self.geometries[Physics.POROUS],
                                        self.time,
                                        self.parameters[str(Physics.POROUS)],
                                        **kwargs)


    def _setup_solvers(self, **kwargs):
        """Set up the solvers for the problems.  
        """

        interval = (self.parameters['start_time'], self.parameters['end_time'])

        if Physics.ELECTRO in self.physics:

            electrosolver_parameters = self.parameters["ElectroSolver"]
            self.electro_solver = SplittingSolver(
                                    self.electro_problem.cardiac_model, 
                                    electrosolver_parameters)

            (vs_, _, _) = self.electro_solver.solution_fields()
            vs_.assign(self.electro_problem.cell_model.initial_conditions())

        if Physics.SOLID in self.physics:
            parameters = self.parameters[str(Physics.SOLID)]
            self.solid_solver = SolidSolver(
                                    self.solid_problem._form, self.time,
                                    interval, parameters['dt'], parameters,
                                    **kwargs)

        if Physics.FLUID in self.physics:
            parameters = self.parameters[str(Physics.FLUID)]
            self.fluid_solver = FluidSolver(
                                    self.fluid_problem._form, self.time,
                                    interval, parameters['dt'], parameters,
                                    **kwargs)

        if Physics.POROUS in self.physics:
            parameters = self.parameters[str(Physics.POROUS)]
            self.porous_solver = PorousSolver(
                                    self.porous_problem._form, self.time,
                                    interval, parameters['dt'], parameters,
                                    **kwargs)


    def _setup_geometries(self, geometry, physics):
        """ Sets up a dictionary with the geometries for the different problems.
        """ 
        self.geometries = {}

        if isinstance(geometry, MultiGeometry):
            for phys in physics:
                try:
                    self.geometries[phys] = geometry.geometries[phys.value]
                except KeyError:
                    msg = "Could not find a geometry for {} physics in "\
                            "MultiGeometry. Ensure that geometry labels "\
                            "correspond to values in Physics "\
                            "enum.".format(phys.value)
                    raise KeyError(msg)

                assert isinstance(self.geometries[phys], HeartGeometry)
        elif len(physics) == 1:
            self.geometries[physics[0]] = geometry
        else:
            for phys in physics:
                self.geometries[phys] = geometry.copy(deepcopy=True)


    def _check_physics_interactions(self):
        # If multiple physics are defined, check that all are involved in an
        # interaction and that physics involved in an interaction are set up
        if len(self.physics) == 0 and len(self.interactions) > 0:
            msg = "At least one interaction has been set up, but no physics\n"\
                    "Interactions: {}".format(self.interactions)
            raise KeyError(msg)
        if len(self.physics) > 1:
            int_physics = set(np.array([ia.to_list()
                                            for ia in self.interactions]).flat)
            for p in int_physics:
                if p not in self.physics:
                    msg = "Physics {} appears in interaction, but is not set "\
                            "up.".format(p)
                    raise KeyError(msg)
            for p in self.physics:
                if p not in int_physics:
                    msg = "Physcis {} is set up, but does not appear in any "\
                            "interaction.".format(p)
                    raise KeyError(msg)


    def update_parameters(self, physics, parameters):
        if physics == Physics.ELECTRO:
            self.parameters.set_electro_parameters(parameters)
        elif physics == Physics.SOLID:
            self.parameters.set_solid_parameters(parameters)
        elif physics == Physics.FLUID:
            self.parameters.set_fluid_parameters(parameters)
        elif physics == Physics.POROUS:
            self.parameters.set_porous_parameters(parameters)