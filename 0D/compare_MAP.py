# -*- coding: utf-8 -*-

from setup0D import create_module
import numpy as np
import matplotlib.pyplot as plt 

dt = 0.01

'''
# nSR
odes = ['hAM_KSMT_nSR', 'FK_nSR']

for i, ode in enumerate(odes):
    module, forward = create_module(ode)
    model_params = module.init_parameter_values()
    states = module.init_state_values()

    V_index = module.state_indices("V")
    offset_index = module.parameter_indices("stim_offset")
    model_params[offset_index] = 10
    
    
    # Start array for recording results
    V = [states[V_index]]

    t = 0; tstop = 400
    while t <= tstop:
        forward(states, t, dt, model_params)
        V.append(states[V_index])
        t += dt

    # Plot the results
    if i==1:
        V = np.array(V)
        V *= 120
        V += -75
        print (-60+75)/120.
    tsteps = np.linspace(0, tstop, len(V))
    plt.plot(tsteps, V, linewidth=1.5)

    
plt.xlabel(r"Time [ms]", fontsize=20)
plt.ylabel(r"V [mV]", fontsize=20)
#plt.title(r'MAP', fontsize=20)
plt.grid()
plt.legend([u"Koivumäki nSR", u"FK nSR"], prop={'size':20})
plt.savefig('../fig/compare_MAPS_nSR.png')
plt.show()
plt.close()

'''
# cAF
odes = ['hAM_KSMT_cAF', 'FK_cAF']

for i, ode in enumerate(odes):
    module, forward = create_module(ode)
    model_params = module.init_parameter_values()
    states = module.init_state_values()

    V_index = module.state_indices("V")
    offset_index = module.parameter_indices("stim_offset")
    model_params[offset_index] = 10
    
    # Start array for recording results
    V = [states[V_index]]

    t = 0; tstop = 400
    while t <= tstop:
        forward(states, t, dt, model_params)
        V.append(states[V_index])
        t += dt

    # Plot the results
    if i==1:
        V = np.array(V)
        V *= 125
        V += -76
        print (-60+75)/120.
    tsteps = np.linspace(0, tstop, len(V))
    plt.plot(tsteps, V, linewidth=1.5)

    
#plt.xlabel(r"Time [ms]", fontsize=20)
#plt.ylabel(r"V [mV]", fontsize=20)
#plt.title(r'MAP', fontsize=20)
plt.grid()
plt.legend([u"Koivumäki cAF", u"Fenton-Karma cAF"], prop={'size':20})
plt.savefig('../fig/compare_MAPS_cAF.png')
plt.show()
