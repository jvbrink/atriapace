from TissueSquare import *

ode = 'hAM_KSMT_nSR'
D = Constant(0.34)
L = 50
amp = -1410*12
plot = {'range_min':-80.0, 'range_max':40.0, 'mode':'color'}
BCL = 280

solver = TissueSquare(ode, D, L=50, stim_amp=amp, plot_args=plot)
solver.plane_wave(BCL, tstop=220*10, liveplot=True)

