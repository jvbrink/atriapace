# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 

odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF', 'FK_nSR', 'FK_cAF']

results = [np.load("../data/results/dynamic_%s.npy" % ode) for ode in odes]

plt.hold('on')
for result in results:
    BCL, ADP, DI = result
    plt.plot(BCL, ADP, linewidth=2.0)

plt.title("Restitution curve, dynamic pacing protocol", fontsize=20)
plt.axis([250, 1050, 120, 240])
plt.xlabel("BCL [ms]", fontsize=20)
plt.ylabel("APD [ms]", fontsize=20)
plt.legend([u"Koivumäki nSR", u"Koivumäki cAF", "Fenton-Karma nSR", "Fenton-Karma cAF"], 'center right', prop={'size':14})
plt.grid()
plt.savefig("../fig/restcurve_APD_dynamic.png")
plt.show()