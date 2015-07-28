""" 
This code solves the monodomain equations on a 2D tissue square.

The equations are solved using an operator splitting scheme,
solving the ODEs using GOSS and the PDEs using FEniCS. 
The cbcbbeat.GOSSplittingSolver is used to keep the implementation
high-level.
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

class TissueSquare:
    def __init__(self, cellmodel, D, L=50, h=0.5, plot_args={},
                 threshold=0.5, stim_amp=-1410*8, scpath="../data/steadycycles/"):
        # Create Dolfin Mesh
        self.L = L; self.h = h
        self.N = int(L/h)
        self.domain = RectangleMesh(0, 0, L, L, self.N, self.N)

        self.D = D

        # Load the cell model from file
        self.ode = cellmodel
        ode = load_ode(cellmodel)
        cellmodel = dolfin_jit(ode, field_states=['V'],
                                    field_parameters=['stim_amplitude', 'stim_offset'])

        # Create the CardiacModel for the given mesh and cell model
        self.t = Constant(0.0)
        heart = CardiacModel(self.domain, self.t, self.D, None, cellmodel)

        # Create the solver and extract the subsolvers
        self.solver = GOSSplittingSolver(heart, self.GOSSparams())
        self.dolfin_solver = self.solver.ode_solver # Solves spatial PDE
        self.ode_solver = self.dolfin_solver._ode_system_solvers[0] # Solves ODE system
        self.scpath = scpath
        self.threshold = threshold
        self.stim_amp = stim_amp
        self.plot_args = plot_args

    def plane_wave(self, BCL, dt=0.1, tstop=0, liveplot=True):
        # Initialize cell by pacing 0D cell model
        try:
            states = np.load(self.scpath+"%s_BCL%d.npy" % (self.ode, BCL))
        except:
            print "Steady cycle at BCL=%d for ODE model: %s not found." % (BCL, self.ode)
            print "Pacing 0D cell model to find it, this may take a minute."
            states = find_steadycycle(self.ode, BCL, dt, odepath="", scpath=self.scpath)
            print "Steady cycle found, proceeding to do 1D tissue simulation."

        # Load steadystate into 2D model
        (u, um) = self.solver.solution_fields()
        for node in range(self.domain.coordinates().shape[0]):
            self.ode_solver.states(node)[:] = states.copy()
            u.vector()[node] = states[0]

        # Set the stimulus parameter field
        V = VectorFunctionSpace(self.domain, 'CG', 1)
        stimulus = Expression(("near(x[0],0)*amp","offset"), amp=self.stim_amp, offset=0)
        stimulus = interpolate(stimulus, V).vector().array()
        
        # Apply the stimulus
        self.ode_solver.set_field_parameters(stimulus)

        d = 10.0
        px = Point(np.array((self.L/2-d, self.L/2)))
        py = Point(np.array((self.L/2+d, self.L/2)))

        xp = 0; yp = 0;
        t0 = 0; t1 = 0;
        if not tstop: tstop = BCL

        for timestep, (u, vm) in self.solver.solve((0, tstop), dt):
            if liveplot: plot(u, **self.plot_args)

            x = u(px)
            y = u(py)

            if not t0 and x > self.threshold:
                t0 = timestep[0]
            if not t1 and y > self.threshold:
                t1 = timestep[0]
                # Calculate conduction velocity
                print "CV: ", 2*d*self.h/(t1-t0)*1e3
            
    def spiral_wave(self, BCL, dt=0.1, tstop=0, S2time=300, savefig=False):
        # Initialize cell by pacing 0D cell model
        try:
            states = np.load(self.scpath+"%s_BCL%d.npy" % (self.ode, BCL))
        except:
            print "Steady cycle at BCL=%d for ODE model: %s not found." % (BCL, self.ode)
            print "Pacing 0D cell model to find it, this may take a minute."
            states = find_steadycycle(self.ode, BCL, dt, odepath="", scpath=self.scpath)
            print "Steady cycle found, proceeding to do 1D tissue simulation."

        # Load steadystate into 2D model
        (u, um) = self.solver.solution_fields()
        for node in range(self.domain.coordinates().shape[0]):
            self.ode_solver.states(node)[:] = states.copy()
            u.vector()[node] = states[0]

        # Set the stimulus parameter fields       
        V = VectorFunctionSpace(self.domain, 'CG', 1)
        S1 = Expression(("near(x[0],0)*amp","offset"), amp=self.stim_amp, offset=0)
        S2 = Expression(("(x[1] < 0.5*L && x[0] < 0.5*L)*amp", "offset"), L=self.L, offset=S2time+1, amp=self.stim_amp)

        S1 = interpolate(S1, V).vector().array()
        S2 = interpolate(S2, V).vector().array()
        nostim = np.zeros(S1.shape, dtype=np.float_)

        cnt = 0

        # Apply S1 stimulus and simulate up until S2
        self.ode_solver.set_field_parameters(S1)
        for timestep, (u, vm) in self.solver.solve((0, S2time), dt):
            fig = plot(u, **self.plot_args)
            if savefig and cnt % 10 == 0:
                padded_index = '%04d' % (cnt/10)
                fig.write_png('tmp/spiral_%s' % self.ode + padded_index)
            cnt += 1
            
        self.ode_solver.set_field_parameters(S2)
        for timestep, (u, vm) in self.solver.solve((S2time, S2time+10), dt):
            fig = plot(u, **self.plot_args)
            if savefig and cnt % 10 == 0:
                padded_index = '%04d' % (cnt/10)
                fig.write_png('tmp/spiral_%s' % self.ode + padded_index)
            cnt += 1
        
        self.ode_solver.set_field_parameters(nostim)
        for timestep, (u, vm) in self.solver.solve((S2time+10, tstop), dt):
            fig = plot(u, **self.plot_args)
            if savefig and cnt % 10 == 0:
                padded_index = '%04d' % (cnt/10)
                fig.write_png('tmp/spiral_%s' % self.ode + padded_index)
            cnt += 1
                     
    def GOSSparams(self):
        params = GOSSplittingSolver.default_parameters()
        params["pde_solver"] = "monodomain"
        params["MonodomainSolver"]["linear_solver_type"] = "iterative"
        params["MonodomainSolver"]["theta"] = 1.0
        params["ode_solver"]["solver"] = "RL1"
        params["apply_stimulus_current_to_pde"] = False
        return params
