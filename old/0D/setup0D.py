# -*- coding: utf-8 -*-

from gotran.codegeneration import compile_module
from gotran.common.options import parameters
from gotran.model.loadmodel import load_ode

def create_module(ode, path="../ode/"):
    """
    Compile a gotran ode and return the module and 
    corresponding forward method.
    """
    generation = parameters.generation.copy()

    for what_not in ["componentwise_rhs_evaluation",
                     "forward_backward_subst",
                     "linearized_rhs_evaluation",
                     "lu_factorization",
                     "jacobian"]:
        generation.functions[what_not].generate = False

    solver = "rush_larsen"
    generation.solvers[solver].generate = True

    if isinstance(ode, str):
        ode = load_ode(path+ode)
        
    module = compile_module(ode, "c", [], generation)
    forward = getattr(module, "forward_"+solver)

    return module, forward


if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt 
    dt = 0.01

    '''
    ode = 'hAM_KSMT_nSR'

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
    plt.plot((100,332), (-60, -60), 'o-',  linewidth=2.0)
    plt.plot((332,600), (-60, -60), 'o-', linewidth=2.0)
    plt.plot((100,600), (-70, -70), 'k--', linewidth=2.0)
    plt.xlabel(r"Time [ms]", fontsize=16)
    plt.ylabel(r"V [mV]", fontsize=16)
    plt.grid()
    plt.legend([u"Action potential", u"APD", u"DI", u"BCL"], 9)
    plt.savefig('../fig/action_potential.png')
    plt.show()
    '''

    BCL = 1000
    odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF']
    for ode in odes:
        # Create module from ode file
        module, forward = create_module(ode)

        # Fetch model parameters and states
        model_params = module.init_parameter_values()
        states = module.init_state_values()

        states = np.load("../data/steadystates/steadystate_%s_BCL%d.npy" % (ode, BCL))

        # Find index of states and parameters of interest
        V_index = module.state_indices("V")
        offset_index = module.parameter_indices("stim_offset")
        BCL_index = module.parameter_indices("stim_period")

        # Set parameters
        model_params[BCL_index] = BCL
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
        plt.plot(tsteps, V, linewidth=1.5)
    
    plt.xlabel(r"Time [ms]", fontsize=16)
    plt.ylabel(r"V [mV]", fontsize=16)
    plt.grid()
    plt.title("Action potential")
    plt.legend([u"Koivumäki nSR", u"Koivumäki cAF"])
    plt.savefig('../fig/action_nSR_cAF_potential.png')
    plt.show()
    