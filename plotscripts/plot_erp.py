# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 

odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF', 'FK_nSR', 'FK_cAF']

results = [np.load("../data/results/erp_%s.npy" % ode) for ode in odes]

plt.hold('on')
for result in results:
    BCL, ERP  = result
    plt.plot(BCL, ERP, linewidth=2.0)

#plt.title("Restitution curve, S1-S2 pacing protocol", fontsize=16)
plt.axis([50, 1050, 50, 415])
plt.xlabel("BCL [ms]", fontsize=20)
plt.ylabel("ERP [ms]", fontsize=20)
plt.legend([u"Koivumäki nSR", u"Koivumäki cAF", "Fenton-Karma nSR", "Fenton-Karma cAF"], 'upper left', prop={'size':16})
plt.grid()
plt.savefig("../fig/erp.png")
plt.show()
