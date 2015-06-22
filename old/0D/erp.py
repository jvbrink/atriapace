from setup0D import create_module
import numpy as np
import matplotlib.pyplot as plt

def erp(ode, BCL_range, dt, start, threshold=0.8):
    """
    Calculate the effective refractatory period of an ode model for
    a range of BCLs. First the cell is paced for 50 cycles at the 
    given S1 BCL, then a S2 pulse of a lower cycle length is issued.
    If the action potential from the S2 pulse has a peak value of 
    less than the threshold ratio of the S1 peak, the cell is refractatory.
    The ERP is the smallest S2 cycle length where the cell is not refractatory.
    """
    module, forward = create_module(ode)
    model_params = module.init_parameter_values()
    init_states = module.init_state_values()
    
    # Find indices of various states/parameters
    V_index = module.state_indices("V")
    BCL_index = module.parameter_indices("stim_period")

    results = np.zeros((2,len(BCL_range)))
    for i, BCL in enumerate(BCL_range):
        # First pace cell for 50 cycles at S1
        try:
            sspath = "../data/steadystates/"
            S1 = np.load(sspath+"steadystate_%s_BCL%d.npy" % (ode, BCL))
        except:
            print "Pacing for BCL %d." % BCL
            # Load initial condition
            states = init_states
            # Set the BCL
            model_params[BCL_index] = BCL
            # Solve ODE system for given number of beats
            t = 0.; tstop = 50*BCL
            while t <= tstop:
                forward(states, t, dt, model_params)
                t += dt
            # Save the final state as the steady state   
            np.save("../data/steadystates/steadystate_%s_BCL%d" % (ode, BCL), states)
            S1 = states.copy()

        if i == 0:
            S2_range = np.arange(start, BCL, 0.1)
        else:
            S2_range = np.arange(results[1][i-1], BCL, 0.1)
        
        for S2 in S2_range:
            print "Trying %g" % S2 
            states = S1.copy()
            model_params[BCL_index] = S2
            V = []
            t = 0; tstop = 2*S2
            while t <= S2:
                forward(states, t, dt, model_params)
                V.append(states[V_index])
                t += dt

            Vpeak = max(V)
            
            V = []
            while t <= tstop:
                forward(states, t, dt, model_params)
                V.append(states[V_index])
                t += dt
            
            if max(V) > 0.8*Vpeak:
                break

        print "BCL: %g  ERP: %s" % (BCL, S2)
        
        results[0][i] = BCL
        results[1][i] = S2

    np.save("../data/results/erp_%s" % ode, results)

if __name__ == '__main__':
    odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF', 'FK_nSR', 'FK_cAF']
    BCL_range = xrange(300, 1005, 5)
    start = [250, 150, 1, 50]
    dt = 0.01


    for i in range(1, len(odes), 3):
        ode = odes[i]
        s = start[i]
        erp(ode, BCL_range, dt, s)
