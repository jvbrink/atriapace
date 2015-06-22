# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 

odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF', 'FK_nSR', 'FK_cAF']

results = [np.load("../data/results/cv_%s.npy" % ode) for ode in odes]

plt.hold('on')
for result in results:
    BCL = result[0][:]
    CV = (result[3][:]+result[2][:]+result[1][:])/3.
    plt.plot(BCL, CV, linewidth=2.0)

plt.title("Restitution curve, conduction velocity")
plt.axis([250, 1050, 660, 800])
plt.xlabel("BCL [ms]", fontsize=16)
plt.ylabel("CV [ms]", fontsize=16)
plt.legend([u"Koivumäki nSR", u"Koivumäki cAF", "Fenton-Karma nSR", "Fenton-Karma cAF"], 'upper right')
plt.grid()
plt.savefig("../fig/restcurve_CV.png")
plt.show()