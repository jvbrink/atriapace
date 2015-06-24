# -*- coding: utf-8 -*-

from gotran.codegeneration import compile_module
from gotran.common.options import parameters
from gotran.model.loadmodel import load_ode
import numpy as np

def create_module(ode, solver="rush_larsen", monitored=[], path="../ode/"):
    """
    Compile a gotran ode and return the module and 
    corresponding forward method.
    """
    generation = parameters.generation.copy()

    # We don't want to compile unessecery features
    for what_not in ["componentwise_rhs_evaluation",
                     "forward_backward_subst",
                     "linearized_rhs_evaluation",
                     "lu_factorization",
                     "jacobian"]:
        generation.functions[what_not].generate = False

    # Select what numerical scheme we want the solver to use
    generation.solvers[solver].generate = True

    if isinstance(ode, str):
        ode = load_ode(path+ode)
    
    module = compile_module(ode, "C", monitored, generation)
    forward = getattr(module, "forward_"+solver)

    return module, forward

def find_steadycycle(ode, BCL, dt, num_of_beats=50, odepath="../ode/", 
                     scpath="../data/steadycycles/"):
    """
    Find a steady cycles by pacing a 0D cell model at a constant BCL
    for a given number of beats. Save the results at the end of a cycle.
    """
    # Compile the ODE solver module
    if isinstance(ode, str):
        module, forward = create_module(ode, path=odepath)
    elif isinstance(ode, list):
        module, forward, ode = ode
    else:
        module, forward = create_module(ode)

    # Fetch model parameter list and init states
    model_params = module.init_parameter_values()
    states = module.init_state_values()

    # Find indices of states/parameters of interest
    index = {}
    index['V'] = module.state_indices('V')
    index['BCL'] = module.parameter_indices("stim_period")

    # Set the BCL of the module
    model_params[index['BCL']] = BCL

    t = 0.; tstop = num_of_beats*BCL
    while t  <= tstop:
        forward(states, t, dt, model_params)
        t += dt

    # Save the final states as the steady cycle
    np.save(scpath+"%s_BCL%d" % (ode, BCL), states)
    return states

def find_steadycycles(ode, BCL_range, dt, num_of_beats=50, odepath="../ode/", 
                     scpath="../data/steadycycles/"):
    """
    Find stady cycles for a range of BCL values for a given ODE.
    """
    # Compile the ODE solver module
    if isinstance(ode, str):
        module, forward = create_module(ode, path=odepath)
    elif isinstance(ode, list):
        module, forward, ode = ode
    else:
        module, forward = create_module(ode)

    # Fetch model parameter list and init states
    model_params = module.init_parameter_values()
    init_states = module.init_state_values()

    # Find indices of states/parameters of interest
    index = {}
    index['V'] = module.state_indices('V')
    index['BCL'] = module.parameter_indices("stim_period")

    # Find steady cycle for each BCL
    for BCL in BCL_range:
        # Load initial conditions
        states = init_states
        # Set the BCL
        model_params[index['BCL']] = BCL

        t = 0.; tstop = num_of_beats*BCL
        while t <= tstop:
            forward(states, t, dt, model_params)
            t += dt

        # Save the final state as the steady cycle
        np.save(scpath+"%s_BCL%d" % (ode, BCL), states)


if __name__ == '__main__':
    ### Example: Create ODE module and plot action potential
    import numpy as np
    import matplotlib.pyplot as plt 
    dt = 0.01

    ode = 'hAM_KSMT_nSR'
    ode = 'hAM_KSMT_cAF'
    ode = 'FK_nSR'
    ode = 'FK_cAF'

    module, forward = create_module(ode)
    model_params = module.init_parameter_values()
    states = module.init_state_values()

    V_index = module.state_indices("V")
    offset_index = module.parameter_indices("stim_offset")
    BCL_index = module.parameter_indices("stim_period")

    # Set parameters
    model_params[BCL_index] = 500
    model_params[offset_index] = 100

    # Start array for recording results
    V = [states[V_index]]

    t = 0; tstop = 700
    while t <= tstop:
        forward(states, t, dt, model_params)
        V.append(states[V_index])
        t += dt

    # Plot the results
    tsteps = np.linspace(0, tstop, len(V))
    plt.plot(tsteps, V, linewidth=2.0)
    plt.xlabel(r"Time [ms]", fontsize=16)
    plt.ylabel(r"V [mV]", fontsize=16)
    plt.grid()
    plt.savefig('../fig/action_potential_%s.png' % ode)
    plt.show()
    
    '''
    ### Example: Find steady cycle of an ODE model

    dt = 0.01
    BCL = 500
    ode = 'hAM_KSMT_nSR'
    find_steadystate(ode, BCL, dt)
    '''

    '''
    ### Example: Find steady cycle for a range of BCLs for same ODE

    dt = 0.01
    BCL_range = range(1000,295,-5)
    ode = 'hAM_KSMT_nSR'
    find_steadystates(ode, BCL_range, dt, num_of_beats=num_of_beats)
    '''