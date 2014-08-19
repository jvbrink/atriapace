from setup0D import create_module
import numpy as np
import matplotlib.pyplot as plt

def pacing_dynamic(ode, BCL_range, dt, pulse_train_length=30e3, 
                   pacing_memory=True, threshold=-60):
    """
    Calculate a restitution curve using dynamic pacing, take 
    in the relevant ode model, the range of BCLs to be used.
    A pulse train of 30 seconds is the default. The regime steps
    straight from one BCL to the next, meaning it contains pacing memory.
    If pacing memory is turned of, we reload the inital states.
    """

    module, forward = create_module(ode)
    model_params = module.init_parameter_values()
    states = module.init_state_values()

    # Find indices of various states/parameters
    V_index = module.state_indices("V")
    BCL_index = module.parameter_indices("stim_period")

    def pulse_train(BCL, init_states):
        """
        Run a pulse train on the cell with given BCL and length, 
        and measure the APD and DI from the last pulse.
        """
        states = init_states
        model_params[BCL_index] = BCL

        t=0.; 
        if pulse_train_length > BCL*50.:
            tstop = np.ceil(pulse_train_length/BCL)*BCL
        else:
            tstop = BCL*50. 

        while t < tstop - BCL:
            forward(states, t, dt, model_params)
            t += dt

        cross_threshold = []
        Vp = states[V_index]
        while t <= tstop:
            forward(states, t, dt, model_params)
            V = states[V_index]
            if (V-threshold)*(Vp-threshold) < 0:
                cross_threshold.append(t)
            Vp = V
            t += dt
        
        APD = cross_threshold[-1] - cross_threshold[-2]
        DI = tstop - cross_threshold[-1] 

        return APD, DI, states

    results_BCL = []
    results_APD = []
    results_DI = []

    prev_states = init_states
    for BCL in BCL_range:
        states = prev_states if pacing_memory else init_states
        APD, DI, prev_states = pulse_train(BCL, states)
        
        results_BCL.append(BCL)
        results_APD.append(APD)
        results_DI.append(DI)

        print "BCL: %g,  APD: %g,   DI: %g" % (BCL, APD, DI)

    np.save("../data/resultss_dynamic_%s" % ode, 
            np.array((results_BCL, results_APD, results_DI)))
    
if __name__ == '__main__':    
    odes = ['FK_nSR', 'FK_cAF']
    thresholds = [0.127720951569, 0.152874943753]

    BCL_range = xrange(1000, 200, -5)
    dt = 0.01

    for i in range(len(odes)):
        ode = odes[i]
        threshold = thresholds[i]
        pacing_dynamic(ode, BCL_range, dt, threshold=threshold)
    
