
from setup0D import create_module, find_steadycycle
import numpy as np
import matplotlib.pyplot as plt

def alterans(ode, BCL_range, dt, threshold=0.1, plot=False,
             odepath='../ode/', scpath="../data/steadycycles/"):
    """
    Calculate a restitution curve using dynamic pacing.
    The ODE model is paced for 50 beats at every BCL value,
    and the APD and DI are measured for the final beat.
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
        try:
            states = np.load(scpath+"%s_BCL%d.npy" % (ode, BCL))
        except:
            print "Steady cycle at BCL=%d for ODE model: %s not found." % (BCL, ode)
            print "Pacing 0D cell model to find it, this may take a minute."
            states = find_steadycycle([module, forward, ode], BCL, dt, odepath=odepath, scpath=scpath)
            print "Steady cycle found, proceeding to do dynamic pacing."

        # Measure when potential crosses threshold
        for i in range(5):
            t = 0; tstop=BCL
            V = [states[index['V']]]
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


if __name__ == '__main__':
    ### Example of use
    ode = 'FK_cAF'
    ode = 'hAM_KSMT_cAF'
    ode = 'FK_nSR'
    ode = 'hAM_KSMT_nSR'

    BCL_range = range(300, 195, -5)
    dt = 0.01

    alterans(ode, BCL_range, dt, plot=True)
