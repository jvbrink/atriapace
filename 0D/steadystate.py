from setup0D import create_module
import numpy as np
import matplotlib.pyplot as plt

def find_steadystate(ode, BCL, dt, num_of_beats=50, path="../ode/"):
    """
    Find a quasi steady state by pacing the 0D cell model at a given cycle length
    for a given number of beats.
    """
    # Compile the ODE solver module
    module, forward = create_module(ode, path=path)

    # Fetch model parameter list and init states
    model_params = module.init_parameter_values()
    states = module.init_state_values()

    # Find indices of various states/parameters
    V_index = module.state_indices("V")
    BCL_index = module.parameter_indices("stim_period")

    # Set the BCL
    model_params[BCL_index] = BCL

    # Start array for recording results
    V = [states[V_index]]

    # Solve ODE system over time
    t = 0.
    tstop = num_of_beats*BCL
    while t <= tstop:
    	forward(states, t, dt, model_params)
    	V.append(states[V_index])
    	t += dt	

    # Save the final state as the steady state
    sspath = "../data/steadystates/"
    np.save(sspath+"steadystate_%s_BCL%d" % (ode, BCL), states)
    return states

def find_steadystates(ode, BCL_range, dt, num_of_beats=50, path="../ode/"):
    """
    Find the quasi steady state for a range of cycle lengths by pacing
    the 0D cell model at a given cycle length for a given number of beats,
    the default being 50 beats.
    """

    # Compile the ODE solver module
    module, forward = create_module(ode, path=path)

    # Fetch model parameter list and init states
    model_params = module.init_parameter_values()
    init_states = module.init_state_values()

    # Find indices of various states/parameters
    V_index = module.state_indices("V")
    BCL_index = module.parameter_indices("stim_period")

    for BCL in BCL_range:
        print "Pacing for BCL %d." % BCL

        # Load initial condition
        states = init_states
        # Set the BCL
        model_params[BCL_index] = BCL

        # Solve ODE system for given number of beats
        t = 0.; tstop = num_of_beats*BCL
        while t <= tstop:
            forward(states, t, dt, model_params)
            t += dt

        # Save the final state as the steady state   
        np.save("../data/steadystates/steadystate_%s_BCL%d" % (ode, BCL), states)

if __name__ == '__main__':
    BCL_range = range(1000,295,-5)
    num_of_beats = 50
    dt = 0.01
    ode = "FK_cAF"

    find_steadystates(ode, BCL_range, dt, num_of_beats=num_of_beats)
