import numpy as np

from dolfin import Constant

from geometry import HeartGeometry, MultiGeometry

from m3h3.setup_parameters import Physics
from m3h3.pde import ElectroProblem, SolidProblem, PorousProblem, FluidProblem
from m3h3.pde.solver import *

from cbcbeat import Expression, Constant

""" 
TODO:
Fix readthedocs, make the code look better.
Fix init files. 
Update initial conditions.
Explain how to add more cell models (lower/upper case)
Should everything be imported from cbcbeat?
Describe the parameters in detail, and where t hey are sent. 
Make it work with 3D? 
What is vs_? Since it can be assigned initial conditions? Is it a bunch of state variables? 
dt should be in solver parameters? 
Rename parameter class to avoid confusion with dolfin imports?
Credit where the cell models come from, or add them straight from cbcbeat instead
of imports? 
Credit cbcbeat for the demos. 
Install vedo for plotting? 
Might also need to have a c++ compiler installed to run meshes from file. 
The read mesh demo should probably be changed a bit? 
Say something about what problem are actually solved. 
Which github link to use?
Make electro parameters class work for M_i and M_e

M3H3 is a framework used for modelling and simulating the heart. It inherits
much of the functionality from other libaries and combines it into a full
framework that encompasses the different branches of physics relevant to
caridac modelling. The goal is for the framework to be able to simulate
the cardiac mechanics, cardiac electrophysiology, hemodynamic modelling, etc.

For now, only methods for simluating the electrical activity is implemented.
For modelling the electrical activity in the heart, the splittingsolver from
cbcbeat is used. Most of the functionality is inherited from the cbcbeat
package. It solves the coupled heart equations, presented in Sundnes et. al.:
https://www.researchgate.net/publication/265487224_Computing_the_Electrical_Activity_in_the_Human_Heart#:~:text=The%20contraction%20of%20the%20heart,called%20the%20electrocardiogram%20(ECG).

"""


class M3H3(object):
    """ A class for representing the full cardiac system.

    This is the class for representing the full cardiac system, and is used to
    do the simulations of the system.

    For the electro problem is solves the Bidomain equations using the
    Splitting solver technique presented in Sundnes (2006). The bidomain
    equations are given by:


    To add stimulus, add it as a keyword argument with the keyword "stimulus".
    If stimulus is dependent on time, use keyword "t" or "time" to connect it
    to the internal timer in the m3h3-object. The stimulus should either be of
    type Expression or Markewise. Use Markerwise if the position of thes
    stimulus should be used as well as if multiple stimuluses should be
    applied.

    The idea is to set up a M3H3 object that can be used to run simulations.
    The basic usage is:

    #. Set up the geometry with a mesh using the fenics.geometry package.
    #. Set the parameters for the cardiac model.
    #. Set the solver parameters.
    #. Create the M3H3 object given the geometry and parameters
    #. Use M3H3 object's step/solve functions for running the simulations.
    #. Do post-processing of the output.

    *Arguments*
        geometry (:py:class:`geometry.Geometry`)
            A geometry object that contains the mesh of the simulation domain.
        parameters (:py:class:`dolfin.Parameters`)
            A Parameters object that contains parameters for setting up the
            cardiac simulation.

    """
    def __init__(self, geometry, parameters, *args, **kwargs):
        self.parameters = parameters
        self.physics = [Physics(p) for p in parameters.keys()
                        if Physics.has_value(p)]
        self.interactions = kwargs.get('interactions', [])
        self.t_step_lengths = self._get_physics_dt()

        if len(self.interactions) > 0:
            self._check_physics_interactions()

        if isinstance(self.parameters["start_time"], Constant) and isinstance(self.parameters["end_time"], Constant):
            self.time = self.parameters['start_time']
            self.parameters["start_time"] = float(self.parameters["start_time"])
            self.parameters["end_time"] = float(self.parameters["end_time"])
        else:
            raise TypeError("start_time and end_time are not df.Constant's")

        self._setup_geometries(geometry, self.physics)
        self._setup_problems(**kwargs)
        self._setup_solvers(**kwargs)

    def step(self):
        """ Does one step for the solvers for a time step equal to the largest
        time step for all the solvers. Finds the interval using the
        internal time variable and the get_max_dt function. Each problem have their own 
        separate time steps. The time steps are assumed to be multiples of each other. 
        This means that if one solver have twice the time step of another solver, 
        the solver with the smallest time step will step two times in the 
        same time interval as the one with the largest time step. This way the 
        solvers will line up at certain times. 

        It then extracts the solvers solutions fields, and updates the m3h3
        objects internal solution fields with it. 

        At the end, it updates m3h3's internal time-variable. 

        """
        # Setup the number of steps for each solver if first iteration.
        if self.time.values()[0] == self.parameters["start_time"]:
            self.num_step, self.max_dt = self._get_num_steps()

        time = float(self.time)

        dt = self._get_physics_dt()

        if Physics.ELECTRO in self.physics:
            for _ in range(self.num_step[str(Physics.ELECTRO)]):
                # Interval to solve for:
                interval = (time, time+dt[str(Physics.ELECTRO)])
                print(interval)

                # Does one step and extracts the solution fields from solver.
                self.electro_solver.step(interval)
                solution_fields = self.electro_solver.solution_fields()
                solution_fields[0].assign(solution_fields[1])

                # Updates m3h3's internal solution fields.
                self.electro_problem.update_solution_fields(solution_fields[0],
                                                            solution_fields[1],
                                                            solution_fields[2])

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

            *Returns*

        """
        self.max_dt = self._get_num_steps()[1]
        t0 = self.parameters["start_time"]
        t1 = t0 + self.max_dt

        while float(self.time) <= self.parameters["end_time"]:
            self.step()
            yield (t0, t1), self.get_solution_fields()
            t0 = t1
            t1 += self.max_dt

    def get_solution_fields(self):
        """ Returns the solution fields for the different problems. Returns a
        dictionary with keys given by strings representing which problem they
        correspond to.

        *Returns*
            solution_fields (:py:class:`dictionary` of :py:class:`dolfin.Function`)
                Returns a dictionary with the solution fields for each problem.
                Keys are the enum value for the problems.  

        *Usage*
            .. code-block:: python 

                sf = system.get_solution_fields()
                sf_electro = sf[Physics.ELECETRO.value]
            
            where system is an instance of m3h3 and the solution fields for 
            the electro problem is returned. 


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
        """ Function for returning the number of steps for each solver

        This function returns the number of step that each solver do. Since 
        m3h3 allows different time-steps for each solver, some solver might 
        have to do multiple steps compared to the other solvers. 

        *Returns*
            num_steps (:py:class:`dictionary` of :py:class:`float`)
                Dictionary of the number of steps for each solver. 
            max_dt (:py:class:`float`)
                The longest time-step for all solvers. 
        """
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
        solvers are multiples of each other. This is important so that they
        match up on certain times. The solvers with shorter time steps are then
        run multiple times between each step of the solvers with longer time 
        steps. See the step function for more information.  

        *Arguments*
            dt (:py:class:`float`)
                Time step
            min_dt (:py:class:`float`)
                Time step
                
        *Return*
            (:py:class:`boolean`)
                Returns True if the input are multiples of each other. 
        
        *Raises*
            Raises an ValueError if the input is not multiples of each other. 

        """
        if not (dt/min_dt).is_integer():
            msg = "Time step sizes have to be multiples of each other."\
                    "{} is not a multiple of {}".format(dt, min_dt)
            raise ValueError(msg)
        else:
            return True

    def _get_physics_dt(self):
        """ Returns a dictionary containing the time steps for each
        solver.
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
                                        **kwargs)

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
        """Function for setting up the solvers.

        This function sets up the solver for each problem with the given parameters
        and keyword arguments. The solver parameters should be in the parameter 
        object. 
        """

        interval = (self.parameters['start_time'], self.parameters['end_time'])

        if Physics.ELECTRO in self.physics:

            electrosolver_parameters = self.parameters["ElectroSolver"]
            self.electro_solver = ElectroSolver(
                                    model = self.electro_problem.cardiac_model,
                                    params = electrosolver_parameters)

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
