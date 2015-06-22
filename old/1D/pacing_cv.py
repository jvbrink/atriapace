"""
Create resitution curves showing conduction velocity against BCL.
Uses gotran ODE models and goss to solve them. The spatial diffusion
model is solved using dolfin and cbcbeat. Uses operator splitting in 
a 1D cable domain.

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

# Import function for pacing 0D cell model to find quasi steady-state
import sys
sys.path.append("../0D/")
from steadystate import find_steadystate

# Turn of dolfin printing
set_log_level(ERROR)

def GOSSparams():
    params = GOSSplittingSolver.default_parameters()
    params["pde_solver"] = "monodomain"
    params["MonodomainSolver"]["linear_solver_type"] = "iterative"
    params["MonodomainSolver"]["theta"] = 1.0
    params["ode_solver"]["solver"] = "RL1"
    params["apply_stimulus_current_to_pde"] = False
    return params

def pacing_cv(ode, BCL_range, D, dt, threshold=0, stim_amp=0, plot_args=None):
    D = Constant(D) # Must be adjusted with temporal and spatial resolution
    L = 20
    h = 0.5
    dt_low = 0.1

    # Create dolfin mesh
    N = int(L/h)
    domain = IntervalMesh(N, 0, L)

    # Load the cell model from file
    ode = load_ode(ode)
    cellmodel = dolfin_jit(ode, field_states=["V"], 
                            field_parameters=["stim_period", "stim_amplitude"])

    # Create the CardiacModel for the given mesh and cell model
    heart = CardiacModel(domain, Constant(0.0), D, None, cellmodel)

    # Create the solver
    solver = GOSSplittingSolver(heart, GOSSparams())

    # Get the solution fields and subsolvers
    dolfin_solver = solver.ode_solver
    ode_solver = dolfin_solver._ode_system_solvers[0]

    results = np.zeros((4, len(BCL_range)))

    for i, BCL in enumerate(BCL_range):
        # Set the stimulus parameter field
        stim_field = np.zeros(2*(N+1), dtype=np.float_)
        stim_field[0] = BCL
        stim_field[1] = stim_amp
        ode_solver.set_field_parameters(stim_field)

        # Pace 0D cell model and find quasi steady state 
        sspath = "../data/steadystates/"
        try:
            states = np.load(sspath+"steadystate_%s_BCL%d.npy" % (ode, BCL))
        except:
            print "Did not find steadystate for %s at BCL: %g, pacing 0D model" % (ode, BCL)
            find_steadystate(ODE, BCL, 0.01)
        
        # Load quasi steady state into 1D model
        (u, um) = solver.solution_fields()
        for node in range(N+1):
            ode_solver.states(node)[:] = states
            u.vector()[node] = states[0] 

        # Used for measuring cell activation
        xp = 0; yp=0
        results[0][i] = BCL
        print "BCL: %g" % BCL
        # Do 3 pulses

        for pulsenr in range(1,4):
            # Solve the pulse in higher temporal resolution
            for timestep, (u, vm) in solver.solve((i*BCL, i*BCL+50), dt):
                if plot_args: plot(u, **plot_args)
                print "Testing"

                x = u.vector()[20]
                y = u.vector()[35]

                if x > threshold and (x-threshold)*(xp-threshold) < 0:
                    t0 = (threshold - xp)/(x - xp)*dt + tp 
                if y > threshold and (y-threshold)*(yp-threshold) < 0:
                    t1 = (threshold - yp)/(y - yp)*dt + tp 
                    # Calculate conduction velocity
                    Cv = 15*h/(t1-t0)*1e3
                    print "\tCv: %g" % Cv
                    results[pulsenr][i] = Cv
                    t0 = 0; t1 = 0

                xp = x
                yp = y
                tp = timestep[0]

            print "adada"

            # Wait for next pulse
            for timestep, (u, vm) in solver.solve((i*BCL+50, (i+1)*BCL), dt_low):
                if plot_args: plot(u, **plot_args)


    np.save("../data/results/cv_%s"%ode, results) 

if __name__ == '__main__':
    odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF', 'FK_nSR', 'FK_cAF']

    thresholds = [0, 0, 0.5, 0.5]
    stim_amps = [-1410*8, -1410*8, -0.8, -0.8]
    Ds = [0.31, 0.31, 0.078, 0.078]

    dt = 0.01
    BCL_range = xrange(1000, 295, -5)

    for i in range(len(odes)-1):
        ode = odes[i]
        threshold = thresholds[i]
        stim_amp = stim_amps[i]
        D = Ds[i]

        plot_args = {'range_min':-80.0, 'range_max':40.0}

        pacing_cv(ode, BCL_range, D, dt, threshold=threshold, stim_amp=stim_amp, plot_args={'range_min':-80., 'range_max':40., 'interactive':True})
    
