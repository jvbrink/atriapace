
from setup0D import create_module, find_steadycycle
import numpy as np
import matplotlib.pyplot as plt

def pacing_S1S2(ode, S1, S2_range, dt, threshold=0.1, offset=10, plot=False,
                odepath='../ode/', scpath="../data/steadycycles/"):
    """
    Calculate a restitution curve using S1-S2 pacing.
    The ODE model is paced for 50 beats at every BCL value,
    and the APD and DI are measured for the final beat.
    """

    # Compile the ODE solver module
    if isinstance(ode, str):
        module, forward = create_module(ode, path=odepath)
    else:
        module, forward, ode = ode

    # Load in inital state after S1 pacing
    try:
        init_states = np.load(scpath+"%s_BCL%d.npy" % (ode, S1))
    except:
        print "Steady cycle at S1=%d for ODE model: %s not found." % (S1, ode)
        print "Pacing 0D cell model to find it, this may take a minute."
        states = find_steadycycle([module, forward, ode], BCL, dt, odepath=odepath, scpath=scpath)
        print "Steady cycle found, proceeding to do S1-S2 pacing."

    # Get state and parameter indices
    model_params = module.init_parameter_values()
    index = {}
    index['V'] = module.state_indices('V')
    index['BCL'] = module.parameter_indices('stim_period')
    
    def pulse(S2):
        """Pulse the steady cycle with a single S2 pulse and measure APD"""
        states = init_states.copy()
        model_params[index['BCL']] = S2

        t = 0
        while t <= S2:
            forward(states, t, dt, model_params)
            t += dt

        V = [states[index['V']]]
        while t <= 2*S2:
            forward(states, t, dt, model_params)
            t += dt
            V.append(states[index['V']])

        # Extract times from results
        V = np.array(V)
        vmin = min(V); vmax = max(V)
        Vthresh = (max(V)-min(V))*threshold + min(V)
        APD = len(V[V > Vthresh])*dt
        DI = S2 - APD
        
        if plot:
            # Plot action potential
            tarray = np.linspace(0, S2, len(V))
            plt.plot(tarray, V, linewidth=1.5)

            # Find intersections and plot them
            above = tarray[V > Vthresh]
            lt, ht = above[0], above[-1]

            plt.plot([lt, ht], [Vthresh, Vthresh], 'o-', linewidth=1.5)
            plt.axis([0, S2+offset, vmin*1.05, vmax*1.05])
            plt.grid()
            plt.xlabel('Time [ms]')
            plt.ylabel('V [rel.]')
            plt.show()

        return APD, DI

    for S2 in S2_range:
        APD, DI = pulse(S2)
        print "S2: %g,\t APD: %g,\t DI: %g" % (S2, APD, DI)

if __name__ == '__main__':
    ### Example of use
    ode = 'hAM_KSMT_nSR'
    ode = 'hAM_KSMT_cAF'
    ode = 'FK_nSR'
    ode = 'FK_cAF'

    dt = 0.01
    S1 = 1000
    S2_range = range(500, 295, -5)

    pacing_S1S2(ode, S1, S2_range, dt, plot=True)
