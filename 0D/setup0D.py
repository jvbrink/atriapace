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

def find_steadycycle(ode, BCL, dt, pace_time=50000, odepath="../ode/", 
                     scpath="../data/steadycycles/"):
    """
    Find a steady cycles by pacing a 0D cell model at a constant BCL
    for a given number of beats. Save the results at the end of a cycle.
    """
    # Compile the ODE solver module
    if isinstance(ode, str):
        module, forward = create_module(ode, path=odepath)
    else:
        module, forward, ode = ode

    # Fetch model parameter list and init states
    model_params = module.init_parameter_values()
    states = module.init_state_values()

    # Find indices of states/parameters of interest
    index = {}
    index['V'] = module.state_indices('V')
    index['BCL'] = module.parameter_indices("stim_period")

    # Set the BCL of the module
    model_params[index['BCL']] = BCL

    num_of_beats = int(np.ceil(pace_time/BCL))
    t = 0.; tstop = num_of_beats*BCL
    while t  <= tstop:
        forward(states, t, dt, model_params)
        t += dt

    # Save the final states as the steady cycle
    np.save(scpath+"%s_BCL%d" % (ode, BCL), states)
    return states

def find_steadycycles(ode, BCL_range, dt, pace_time=50000, odepath="../ode/", 
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
        print "Pacing %s at BCL: %g..." % (ode, BCL),
        # Load initial conditions
        states = init_states
        # Set the BCL
        model_params[index['BCL']] = BCL

        num_of_beats = int(np.ceil(pace_time/BCL))
        t = 0.; tstop = num_of_beats*BCL
        while t <= tstop:
            forward(states, t, dt, model_params)
            t += dt

        # Save the final state as the steady cycle
        np.save(scpath+"%s_BCL%d" % (ode, BCL), states)
        print "Done"


if __name__ == '__main__':
    ### Example: Find steady cycle for a range of BCLs for same ODE

    dt = 0.01
    BCL_range = range(1000, 295, -5)
    ode = 'hAM_KSMT_cAF'
    find_steadycycles(ode, BCL_range, dt)
