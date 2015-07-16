# -*- coding: utf-8 -*-

from setup0D import create_module, find_steadycycle
import numpy as np
import matplotlib.pyplot as plt 

def plot_ap(ode, BCL=1000, dt=0.1, offset=10, savefig='', odepath="../ode/"):
    for i, ode in enumerate(odes):
        # Compile the ODE solver module
        module, forward = create_module(ode, path=odepath)
        
        # Fetch model parameter list and init states
        model_params = module.init_parameter_values()
        
        # Find indices of states/parameters of interest
        index = {}
        index['V'] = module.state_indices('V')
        index['BCL'] = module.parameter_indices("stim_period")
        index['offset'] = module.parameter_indices("stim_offset")

        # Load in inital state from steadycycle
        try:
            states = np.load(scpath+"%s_BCL%d.npy" % (ode, BCL))
        except:
            print "Steady cycle at BCL=%d for ODE model: %s not found." % (BCL, ode)
            print "Pacing 0D cell model to find it, this may take a minute."
            states = find_steadycycle([module, forward, ode], BCL, dt, odepath=odepath, scpath=scpath)
            print "Steady cycle found, proceeding to simulate action potential."

        # Set ODE parameters
        model_params[index['BCL']] = BCL
        model_params[index['offset']] = offset

        # Start array for recording results
        V = [states[index['V']]]

        t = 0; tstop = BCL+offset
        while t <= tstop:
            forward(states, t, dt, model_params)
            V.append(states[index['V']])
            t += dt

        # Rescale the results
        # V = np.array(V)
        # V -= min(V)
        # V /= np.max(V)

        # Plot the results
        t = np.linspace(0, tstop, len(V))
        plt.plot(t, V, linewidth= 1.5)

def compare_ap(odes, BCL, dt, offset=10, savefig='', odepath='../ode/', scpath="../data/steadycycles/"):
    '''
    Compare the action potential of two different ODE cell models at a 
    steady cycle with a given BCL.
    '''
    for i, ode in enumerate(odes):
        # Compile the ODE solver module
        module, forward = create_module(ode, path=odepath)
        
        # Fetch model parameter list and init states
        model_params = module.init_parameter_values()
        
        # Find indices of states/parameters of interest
        index = {}
        index['V'] = module.state_indices('V')
        index['BCL'] = module.parameter_indices("stim_period")
        index['offset'] = module.parameter_indices("stim_offset")

        # Load in inital state from steadycycle
        try:
            states = np.load(scpath+"%s_BCL%d.npy" % (ode, BCL))
        except:
            print "Steady cycle at BCL=%d for ODE model: %s not found." % (BCL, ode)
            print "Pacing 0D cell model to find it, this may take a minute."
            states = find_steadycycle([module, forward, ode], BCL, dt, odepath=odepath, scpath=scpath)
            print "Steady cycle found, proceeding to simulate action potential."

        # Set ODE parameters
        model_params[index['BCL']] = BCL
        model_params[index['offset']] = offset

        # Start array for recording results
        V = [states[index['V']]]

        t = 0; tstop = BCL+offset
        while t <= tstop:
            forward(states, t, dt, model_params)
            V.append(states[index['V']])
            t += dt

        # Rescale the results
        # V = np.array(V)
        # V -= min(V)
        # V /= np.max(V)

        # Plot the results
        t = np.linspace(0, tstop, len(V))
        plt.plot(t, V, linewidth= 1.5)

    plt.xlabel(r"Time [ms]", fontsize=20)
    plt.ylabel(r"V [rel.]", fontsize=20)
    plt.title(r'MAP', fontsize=20)
    plt.axis([0, tstop, -0.1, 1.1])
    plt.grid()
    plt.legend(odes, prop={'size':20})
    plt.show()
    if savefig: plt.savefig(savefig)
    plt.close()

if __name__ == '__main__':
    odes_nSR = ['FK_nSR', 'broken']
    odes_cAF = ['hAM_KSMT_cAF', 'FK_cAF']
    BCL = 500
    dt = 0.01

    compare_ap(odes_nSR, BCL, dt)
#    compare_ap(odes_cAF, BCL, dt)

