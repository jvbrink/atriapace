"""
This code solves the monodomain equation for a given ode system on a 1D 
tissue strand using cbcbeat and operator splitting.

Note on stimulus: The stimulus can be controlled through the cbcbeat
system and solver. In this implementation however, it is controlled directly
through the cellmodel parameters (i.e., through the .ode file). 

Protocol:
First a 0D cell model is paced for 50 cycles to reach a quasi steady-state,
this state is loaded into a 1D tissue strand, and 5 pulses are initiated
through the strand. For the last pulse, the conduction velocity is measured
The diffusion constant D has been adjusted so that the nSR model has a velocity
of roughly 750 mm/s at a BCL of 1000 ms.
"""

from dolfin import *
from gotran import load_ode
from goss import dolfin_jit
from cbcbeat.cardiacmodels import CardiacModel
from cbcbeat.gossplittingsolver import GOSSplittingSolver
import numpy as np
import matplotlib.pyplot as plt

# Import function for pacing 0D cell model to find quasi steady-state
import sys
sys.path.append("../0D/")
from setup0D import create_module, find_steadycycle

# Turn of dolfin printing
set_log_level(ERROR)

class TissueStrand:
    def __init__(self, cellmodel, D, L=20, h=0.5, plot_args={},
                 threshold=0, stim_amp=-1410*8, scpath="../data/steadycycles/"):
        self.D = Constant(D)
        self.L = 20
        self.h = 0.5

        # Create dolfin mesh
        self.N = int(L/h)
        domain = IntervalMesh(self.N, 0, self.L)

        # Load the cell model from file
        self.ode = cellmodel
        ode = load_ode(cellmodel)
        cellmodel = dolfin_jit(ode, field_states=['V'],
                                    field_parameters=['stim_period', 'stim_amplitude'])

        # Create the CardiacModel for the given mesh and cell model
        t = Constant(0.0)
        heart = CardiacModel(domain, t, self.D, None, cellmodel)

        # Create the solver and extract the subsolvers
        self.solver = GOSSplittingSolver(heart, self.GOSSparams())
        self.dolfin_solver = self.solver.ode_solver # Solves spatial PDE
        self.ode_solver = self.dolfin_solver._ode_system_solvers[0] # Solves ODE system
        self.scpath = scpath
        self.threshold = threshold
        self.stim_amp = stim_amp
        self.plot_args = plot_args

    def GOSSparams(self):
        params = GOSSplittingSolver.default_parameters()
        params["pde_solver"] = "monodomain"
        params["MonodomainSolver"]["linear_solver_type"] = "iterative"
        params["MonodomainSolver"]["theta"] = 1.0
        params["ode_solver"]["solver"] = "RL1"
        params["apply_stimulus_current_to_pde"] = False
        return params

    def pulse(self, BCL, dt=0.01, num_of_pulses=5, threshold=0, liveplot=False):
        # Initialize cell by pacing 0D cell model
        try:
            states = np.load(self.scpath+"%s_BCL%d.npy" % (self.ode, BCL))
        except:
            print "Steady cycle at BCL=%d for ODE model: %s not found." % (BCL, self.ode)
            print "Pacing 0D cell model to find it, this may take a minute."
            states = find_steadycycle(self.ode, BCL, dt, scpath=self.scpath)
            print "Steady cycle found, proceeding to do 1D tissue simulation."

        # Set the stimulus parameter field
        N, L, h = self.N, self.L, self.h
        stim_field = np.zeros(2*(N+1), dtype=np.float_)
        stim_field[::2] = BCL
        stim_field[-1] = self.stim_amp
        self.ode_solver.set_field_parameters(stim_field)

        # Load states into 1D model
        (u, um) = self.solver.solution_fields()
        for node in range(N+1):
            self.ode_solver.states(node)[:] = states.copy()
            u.vector()[node] = states[0]

        # Do 5 pulses and measure CV
        threshold = self.threshold
        for i in range(5):
            t0 = 0; t1 = 0
            for timestep, (u, vm) in self.solver.solve((i*BCL, i*BCL+50), dt):
                if liveplot: plot(u, **self.plot_args)

                x = u.vector()[35]; y = u.vector()[20]
                if x > threshold and (x-threshold)*(xp-threshold) < 0:
                    t0 = (threshold - xp)/(x - xp)*dt + tp 
                if y > threshold and (y-threshold)*(yp-threshold) < 0:
                    t1 = (threshold - yp)/(y - yp)*dt + tp 
                    Cv = 15*h/(t1-t0)*1e3
                    print "BCL: %d\t Pulse: %d\t CV: %g" % (BCL, i+1, Cv)
                    t0 = 0; t1 = 0

                xp = x
                yp = y
                tp = timestep[0]

            for timestep, (u, vm) in self.solver.solve((i*BCL+50, (i+1)*BCL), 10*dt):
                if liveplot: plot(u, **self.plot_args)

if __name__ == '__main__':
    # Only one of these should be 'active' at once
    #solver = TissueStrand('hAM_KSMT_nSR', 0.31, stim_amp=-1410*8, threshold=0, plot_args={'range_min':-80.0, 'range_max':40.0})
    #solver = TissueStrand('hAM_KSMT_cAF', 0.31, stim_amp=-1410*8, threshold=0, plot_args={'range_min':-80.0, 'range_max':40.0})
    solver = TissueStrand('FK_nSR', 0.077, stim_amp=-0.8, threshold=0.5, plot_args={'range_min':0.0, 'range_max':1.0})
    #solver = TissueStrand('FK_cAF', 0.077, stim_amp=-0.8, threshold=0.5, plot_args={'range_min':0.0, 'range_max':1.0})

    solver.pulse(1000, num_of_pulses=5, liveplot=True)

