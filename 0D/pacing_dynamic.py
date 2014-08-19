from setup0D import create_module
import numpy as np
import matplotlib.pyplot as plt

def pacing_dynamic(ode, BCL_range, dt, threshold=-60):
    """
    Calculate a restituion curve using dynamic pacing. The ODE model 
    is paced for 50 beats at every BCL value, and APD and DI are measured
    for the final beat.
    """

    # Create ODE solver
    module, forward = create_module(ode)

    # Get model parameters and initial conditions
    model_params = module.init_parameter_values()
    init_states = module.init_state_values()

    # Get parameter indices
    V_index = module.state_indices("V")
    BCL_index = module.parameter_indices("stim_period")

    # For storing the results
    results = np.zeros((3, len(BCL_range)))

    for i, BCL in enumerate(BCL_range):
        # Set BCL
        model_params[BCL_index] = BCL

        # Read in steady state from file
        try:
            sspath = "../data/steadystates/"
            states = np.load(sspath+"steadystate_%s_BCL%d.npy" % (ode, BCL))
        except:
            print "Pacing for BCL %d." % BCL
            # Load initial condition
            states = init_states
            # Solve ODE system for given number of beats
            num_of_beats = 50
            t = 0.; tstop = num_of_beats*BCL
            while t <= tstop:
                forward(states, t, dt, model_params)
                t += dt
            # Save the final state as the steady state   
            np.save(sspath+"steadystate_%s_BCL%d" % (ode, BCL), states)

        # Used to measure when potential crosses threshold
        cross_threshold = []
        Vp = states[V_index]

        t=0.; tstop = BCL;
        while t < tstop:
            forward(states, t, dt, model_params)
            V = states[V_index]

            if (V-threshold)*(Vp-threshold) < 0:
                cross_threshold.append(t)
            Vp = V
            t += dt

        try:
            APD = cross_threshold[-1] - cross_threshold[-2]
            DI = tstop - cross_threshold[-1]
        except:
            APD = 0
            DI = 0

        print "BCL: %g,  APD: %g,   DI: %g" % (BCL, APD, DI)

        results[0][i] = BCL
        results[1][i] = APD
        results[2][i] = DI

    np.save("../data/results/dynamic_%s" % ode, results)

if __name__ == '__main__':    
    odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF', 'FK_nSR', 'FK_cAF']
    thresholds = [-60, -60, 0.127720951569, 0.152874943753]

    BCL_range = xrange(1000, 295, -5)
    dt = 0.01

    for i in range(1,len(odes)):
        ode = odes[i]
        threshold = thresholds[i]
        pacing_dynamic(ode, BCL_range, dt, threshold=threshold)
    
