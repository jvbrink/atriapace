"""
This script uses the TissueSquare class to initiate a spiral wave 
in a 2D square tissue. The spiral wave plots are saved to the 
/tmp/ directory and can be made into a movie when the simulation
is finished.
"""

from TissueSquare import *
 
# Only one of these should be 'active' at once
solver = TissueSquare('hAM_KSMT_nSR', Constant(0.34), L=100, stim_amp=-1410*8, plot_args={'range_min':-80.0, 'range_max':40.0, 'mode':'color'})
#solver = TissueSquare('hAM_KSMT_cAF', Constant(0.34), L=100, stim_amp=-1410*8, plot_args={'range_min':-80.0, 'range_max':40.0, 'mode':'color'})
#solver = TissueSquare('FK_nSR', Constant(0.34), L=100, stim_amp=-0.8, plot_args={'range_min':0.0, 'range_max':1.0, 'mode':'color'})
#solver = TissueSquare('FK_cAF', Constant(0.34), L=100, stim_amp=-0.8, plot_args={'range_min':0.0, 'range_max':1.0, 'mode':'color'})

solver.spiral_wave(500, S2time=200, tstop=1000, dt=0.1, savefig=True)

