
from setup0D import create_module, find_steadycycle
import numpy as np
import matplotlib.pyplot as plt

def find_APD(ode, BCL_range, dt=0.01, threshold=0.1,
             plot=False, odepath='../ode/', scpath='../data/steadycycles/',
             use_prepaced_data=True):
    """
    Calculate the action potential duration (APD) for a given range of basic
    cycle length (BCL) values. The threshold for APD can be adjusted, default
    is 0.1 (APD90). The cell is paced for 50 seconds before the actual measurement
    takes place.
    """
    # Compile the ODE solver module
    if isinstance(ode, str):
        module, forward = create_module(ode, path=odepath)
    else:
        module, forward, ode = ode

    # Get model parameters and initial conditions
    model_params = module.init_parameter_values()
    init_states = module.init_state_values()

    # Get state and parameter indices
    index = {}
    index['V'] = module.state_indices('V')
    index['BCL'] = module.parameter_indices('stim_period')

    # For storing the results
    results = np.zeros((3, len(BCL_range)))

    for i, BCL in enumerate(BCL_range):
        # Set BCL
        model_params[index['BCL']] = BCL

        # Read in steady cycle from file 
        if use_prepaced_data:
            try:
                states = np.load(scpath+"%s_BCL%d.npy" % (ode, BCL))
            except:
                print "Steady cycle at BCL=%d for ODE model: %s not found." % (BCL, ode)
                print "Pacing 0D cell model to find it, this may take a minute."
                states = find_steadycycle([module, forward, ode], BCL, dt, odepath=odepath, scpath=scpath)
                print "Steady cycle found, proceeding to do dynamic pacing."
        else:
            states = init_states.load()


        # Measure when potential crosses threshold
        t = 0; tstop = BCL
        V = [states[index['V']]]

        # First beat
        while t < tstop:
            forward(states, t, dt, model_params)
            t += dt
            V.append(states[index['V']])

        # Extract times from results
        V = np.array(V)
        vmin = min(V); vmax = max(V)
        Vthresh = (max(V)-min(V))*threshold + min(V)
        APD = len(V[V > Vthresh])*dt
        DI = BCL - APD

        print "BCL: %g,  APD: %g,   DI: %g" % (BCL, APD, DI)

        if plot:
            # Plot action potential
            tarray = np.linspace(0, BCL, len(V))
            plt.plot(tarray, V, linewidth=1.5)

            # Find intersections and plot them
            above = tarray[V > Vthresh]
            lt, ht = above[0], above[-1]

            plt.plot([lt, ht], [Vthresh, Vthresh], 'o-', linewidth=1.5)
            plt.axis([0, BCL, vmin*1.05, vmax*1.05])
            plt.grid()
            plt.xlabel('Time [ms]')
            plt.ylabel('V [rel.]')
            plt.show()

if __name__ == '__main__':
    ### Example of use
    ode = 'FK_cAF'
    ode = 'hAM_KSMT_cAF'
    ode = 'FK_nSR'
    ode = 'hAM_KSMT_nSR'

    BCL_range = range(1000, 295, -5)
    dt = 0.01

    pacing_dynamic(ode, BCL_range, dt, plot=True)
