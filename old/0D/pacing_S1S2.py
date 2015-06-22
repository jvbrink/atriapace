from setup0D import create_module
import numpy as np
import matplotlib.pyplot as plt

# Import function for pacing 0D cell model to find quasi steady-state
import sys
sys.path.append("../0D/")
from steadystate import find_steadystate

def pacing_S1S2(ode, S1, S2_range, dt, threshold=-60):
    """
    Calculate a restitution curve using S1-S2 pacing, take in
    the relevant ode model, the S1 length and the range of S2
    values to calculate APD/DI. Note that the steady state of 
    a BCL=S1 must be lying in /data/. If this does not exist 
    for the given ODE file, you can make it using steadystate.py.
    """
    module, forward = create_module(ode)
    model_params = module.init_parameter_values()
    
    # Pace 0D cell model and find quasi-steady state
    try:
        sspath = "../data/steadystates/"
        init_states = np.load(sspath+"steadystate_%s_BCL%d.npy" % (ode, S1))
    except:
        print "Pacing 0D model for BCL %g..." % S1,
        init_states = find_steadystate(S1, 50, dt, ode, plot_results=False)
        print "done."

    # Find indices of various states/parameters
    V_index = module.state_indices("V")
    BCL_index = module.parameter_indices("stim_period")

    def pulse_S2(S2):
        """Pulse the steady state with a single S2 pulse and measure APD/DI."""
        states = init_states
        model_params[BCL_index] = S2

        t = 0.; tstop = 2*S2
        while t <= S2:
            forward(states, t, dt, model_params)
            t += dt

        Vp = states[V_index]
        cross_threshold = []
        while t <= tstop:
            forward(states, t, dt, model_params)
            V = states[V_index]
            if (Vp-threshold)*(V-threshold) < 0:
                cross_threshold.append(t)
            Vp = V
            t += dt

        try:
            APD = cross_threshold[-1] - cross_threshold[-2]
            DI = tstop - cross_threshold[-1]
        except:
            APD = 0
            DI = 0

        return APD, DI

    results = np.zeros((3, len(S2_range)))

    for S2 in S2_range:
        APD, DI = pulse_S2(S2)
        results[0][i] = S2
        results[1][i] = APD
        results[2][i] = DI

        print "S2: %g,  APD: %g,   DI: %g" % (S2, APD, DI)

    np.save("../data/results/S1S2_%s" % ode, results)

if __name__ == '__main__':
    odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF', 'FK_nSR', 'FK_cAF']
    thresholds = [-60, -60, 0.127720951569, 0.152874943753]

    S1 = 1000
    S2_range = xrange(1000, 210, -5)
    dt = 0.01

    for i in range(1):
        ode = odes[i]
        threshold = thresholds[i]
        pacing_S1S2(ode, S1, S2_range, dt, threshold=threshold)