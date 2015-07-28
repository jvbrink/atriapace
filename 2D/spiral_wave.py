"""
This script uses the TissueSquare class to initiate a spiral wave 
in a 2D square tissue. The spiral wave plots are saved to the 
/tmp/ directory and can be made into a movie when the simulation
is finished.
"""

from TissueSquare import *
 
def run_plane_wave_hAM_KSMT_nSR():
    ode = 'hAM_KSMT_nSR'
    D = Constant(0.34)
    L = 100
    amp = -1410*8
    plot = {'range_min':-80.0, 'range_max':40.0, 'mode':'color'}
    BCL = 500

    solver = TissueSquare(ode, D, L=L, stim_amp=amp, plot_args=plot)
    solver.spiral_wave(BCL, S2time=200, tstop=1000, dt=0.1, savefig=True)

def run_plane_wave_hAM_KSMT_cAF():
    ode = 'hAM_KSMT_cAF'
    D = Constant(0.34)
    L = 100
    amp = -1410*8
    plot = {'range_min':-80.0, 'range_max':40.0, 'mode':'color'}
    BCL = 500

    solver = TissueSquare(ode, D, L=L, stim_amp=amp, plot_args=plot)
    solver.spiral_wave(BCL, S2time=200, tstop=1000, dt=0.1, savefig=True)

def run_plane_wave_FK_nSR():
    ode = 'FK_nSR'
    D = Constant(0.26) # Dependant on both dt and h
    L = 100
    amp = -0.8
    plot = {'range_min':0.0, 'range_max':1.0, 'mode':'color'}
    BCL = 500

    solver = TissueSquare(ode, D, L=L, stim_amp=amp, plot_args=plot)
    solver.spiral_wave(BCL, S2time=200, tstop=1000, dt=0.1, savefig=True)

def run_plane_wave_FK_cAF():
    ode = 'FK_cAF'
    D = Constant(0.26) 
    L = 100
    amp = -0.8
    plot = {'range_min':0.0, 'range_max':1.0, 'mode':'color'}
    BCL = 500

    solver = TissueSquare(ode, D, L=L, stim_amp=amp, plot_args=plot)
    solver.spiral_wave(BCL, S2time=200, tstop=1000, dt=0.1, savefig=True)

run_plane_wave_hAM_KSMT_nSR()
#run_plane_wave_hAM_KSMT_cAF()
#run_plane_wave_FK_nSR()
#run_plane_wave_FK_cAF()