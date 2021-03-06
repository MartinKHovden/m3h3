Demo with a complex stimulus
===================================================== 

.. code-block:: python 

    """ 
    A simple example on how to use the m3h3 module for solving the 
    bidomina/monodomain equations coupled to a cell model with a simple domain 
    and a complex stimulus.  

    This example shows how to:
    - Set up a simple geometry using functionality from fenics. 
    - Update the parameters of the electro problem.
    - Update the solver parameters for the electro solver.
    - Add a more complex stimulus using Meshfunctions and Compiled subdomains. 
    - Run the electro simulation.
    - Plot the results.

    """

    from m3h3 import *
    import matplotlib.pyplot as plt 

    # Set up the mesh:
    mesh = UnitSquareMesh(100, 100)

    # Set up the stimulus subdomains: 
    stimulus_domain = MeshFunction("size_t", mesh, mesh.topology().dim())
    stimulus_domain.set_all(0)

    stimulus_1 = CompiledSubDomain("pow(x[0],2) + pow(x[1],2) <= 0.5 + tol", tol = 1e-15 )
    stimulus_1.mark(stimulus_domain, 1)

    stimulus_2 = CompiledSubDomain("pow(x[0]-0.5,2) + pow(x[1]-1.0, 2) <= 0.1 + tol", tol = 1e-15)
    stimulus_2.mark(stimulus_domain, 2)

    # Set up the geometry given the computational domain: 
    geo = Geometry2D(mesh)

    # Set up dt, t_0, and t_max: 
    dt = 0.1
    start_time = Constant(0.0)
    end_time = Constant(1.0)
    num_steps = int((float(end_time) - float(start_time))/dt)

    # Define the conductivity (tensors):
    M_i = 2.0
    M_e = 1.0

    # Set up the parameteres for the heart-model: 
    params = Parameters("M3H3")

    params["end_time"] =end_time
    params["start_time"] =start_time

    # Set up various parameters for the stimulus:
    duration = 2. # ms
    chi = 140.0     # mm^{-1}
    # Membrane capacitance
    C_m = 0.01 # mu F / mm^2
    A = 50000. # mu A/cm^3
    cm2mm = 10.
    factor = 1.0/(chi*C_m) # NB: cbcbeat convention
    amplitude = factor*A*(1./cm2mm)**3 # mV/ms

    I_s = Expression("t >= start ? (t <= (duration + start) ? amplitude : 0.0) : 0.0",
                    t=start_time,
                    start=0.0,
                    duration=duration,
                    amplitude=amplitude,
                    degree=0)

    # Set up Markerwise object for stimulus: 
    stimulus = Markerwise((I_s, I_s), (1,2), stimulus_domain)

    params.set_electro_parameters()

    electro_params = params["Electro"]
    electro_params["dt"] = dt
    electro_params["M_i"] = M_i
    electro_params["M_e"] = M_e
    electro_params["cell_model"]  = "Beeler_reuter_1977"#"Tentusscher_panfilov_2006_M_cell"
    electro_params["stimulus"]= stimulus
    electro_params["applied_current"] = None

    # Set up the parameters for the splitting solver: 
    electrosolver_parameters = params["ElectroSolver"]
    electrosolver_parameters["theta"] = 0.5                        # Second order splitting scheme
    electrosolver_parameters["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
    electrosolver_parameters["CardiacODESolver"]["scheme"] = "RL1" # 1st order Rush-Larsen for the ODEs
    electrosolver_parameters["MonodomainSolver"]["linear_solver_type"] = "iterative"
    electrosolver_parameters["MonodomainSolver"]["algorithm"] = "cg"
    electrosolver_parameters["MonodomainSolver"]["preconditioner"] = "sor"#"petsc_amg"
    electrosolver_parameters["apply_stimulus_current_to_pde"] = True

    # Initialize the system with parameters and geometry.
    system = M3H3(geo, params)

    # Run the simulation by using the step function:
    for i in range(num_steps):
        print("Time interval: (%.2f, %.2f)" % (float(system.time), float(system.time) + dt) )
        system.step()

    # Extract the solution:
    vs_, vs, vur = system.get_solution_fields()[str(Physics.ELECTRO)]

    # Plot the resulting solution fields:
    plt.figure()
    plot(vs[0], title="Transmembrane potential (v) at end time")
    plt.savefig("TransmembranePot.png")
    plt.figure()
    plot(vs[-1], title="1st state variable (s_0) at end time")
    plt.savefig("s_0(T).png")

    print("Done!!")