from dolfin import *
from gotran import load_ode
from goss import dolfin_jit
from beatadjoint.cardiacmodels import CardiacModel
from beatadjoint.gossplittingsolver import GOSSplittingSolver
import numpy as np

# Turn of dolfin printing
set_log_level(ERROR)

def GOSSparams():
    params = GOSSplittingSolver.default_parameters()
    params["pde_solver"] = "monodomain"
    params["MonodomainSolver"]["linear_solver_type"] = "iterative"
    params["MonodomainSolver"]["theta"] = 1.0
    params["ode_solver"]["scheme"] = "RL1"
    params["apply_stimulus_current_to_pde"] = False
    return params

def spiral(ode, domain, D, S1, S2, plot_args={}, save_fig=False, threshold=0, BCL=500):
    # Load the cell model from file
    cellmodel = load_ode(ode)
    cellmodel = dolfin_jit(cellmodel, field_states, field_parameters)

    # Create the CardiacModel for the given domain and cell model
    heart = CardiacModel(domain, Constant(0.0), D, None, cellmodel)

    # Create the solver
    solver = GOSSplittingSolver(heart, GOSSparams())
    dolfin_solver = solver.ode_solver
    ode_solver = dolfin_solver._ode_system_solvers[0]

    # Calculate parameter fields
    V = VectorFunctionSpace(domain, 'CG', 1)
    S1 = interpolate(S1, V).vector().array()
    S2 = interpolate(S2, V).vector().array()
    nostim = np.zeros(S1.shape, dtype=np.float_)

    # Read in steady state from file
    try:
        sspath = "../data/steadystates/"
        steadystate = np.load(sspath+"steadystate_%s_BCL%d.npy" % (ode, BCL))
    except:
        print "Did not find steadystate for %s at BCL: %g, pacing 0D model" % (ode, BCL)
        steadystate = find_steadystate(ODE, BCL, 0.01)

    # Load steadystate into 2D model
    (u, um) = solver.solution_fields()
    for node in range(domain.coordinates().shape[0]):
        ode_solver.states(node)[:] = steadystate
        u.vector()[node] = steadystate[0]

    cnt = 0
    # Stimulate S1 pulse
    plot(u, interactive=True, **plot_args)
    ode_solver.set_field_parameters(S1)
    for timestep, (u, vm) in solver.solve((0, S2_time), dt):
        fig = plot(u, interactive=False, **plot_args)

        if save_fig and cnt % 5 == 0:
            padded_index = '%08d' % cnt
            fig.write_png('../tmp/spiral_%s' % ode + padded_index)
        cnt += 1

    ode_solver.set_field_parameters(S2)

    # Stimulate S2 pulse
    for timestep, (u, vm) in solver.solve((S2_time, S2_time+100), dt):
        fig = plot(u, interactive=False, **plot_args)

        if save_fig and cnt % 5 == 0:    
            padded_index = '%08d' % cnt
            fig.write_png('../tmp/spiral_%s' % ode + padded_index)
        cnt += 1

    # Turn of all stimulus and iterate until simulation stop
    # We monitor the top right corner to find the spiral period
    ode_solver.set_field_parameters(nostim)
    trig_time = []
    xp = 0;
    for timestep, (u, vm) in solver.solve((S2_time+100, tstop), dt):
        fig = plot(u, interactive=False, **plot_args)

        x = u.vector()[-1]
        if x > threshold and (x-threshold)*(xp-threshold) < 0:
            trig_time.append(timestep[0])
            print "Trigger: %g" % timestep[0]
        xp = x

        if save_fig and cnt % 5 == 0:
            padded_index = '%08d' % cnt
            fig.write_png('../tmp/spiral_%s' % ode + padded_index)
        cnt += 1


def planar_wave(ode, domain, D, S1, threshold, d=5):
    '''
    Finds the conduction velocity of a planar wave. First paces the 0D ode
    model for 50 cycles at a BCL of 1000, then initiates a single planar
    wave, measuring the conduction velocity half-way in the domain.
    '''
    # Load the cell model from file
    cellmodel = load_ode(ode)
    cellmodel = dolfin_jit(cellmodel, field_states, field_parameters)

    # Create the CardiacModel for the given domain and cell model
    heart = CardiacModel(domain, Constant(0.0), D, None, cellmodel)

    # Create the solver
    solver = GOSSplittingSolver(heart, GOSSparams())
    dolfin_solver = solver.ode_solver
    ode_solver = dolfin_solver._ode_system_solvers[0]

    BCL = 1000
    try:
        sspath = "../data/steadystates/"
        steadystate = np.load(sspath+"steadystate_%s_BCL%d.npy" % (ode, BCL))
    except:
        print "Did not find steadystate for %s at BCL: %g, pacing 0D model" % (ode, BCL)
        steadystate = find_steadystate(ODE, BCL, 0.01)

    # Load steadystate into 2D model
    (u, um) = solver.solution_fields()
    for node in range(domain.coordinates().shape[0]):
        ode_solver.states(node)[:] = steadystate
        u.vector()[node] = steadystate[0]
    
    # Define the planar wave stimulus
    V = VectorFunctionSpace(domain, 'CG', 1)
    S1 = interpolate(S1, V).vector().array()

    # Apply the stimulus
    ode_solver.set_field_parameters(S1)

    xp = 0; yp = 0;
    t0 = 0; t1 = 0;
    for timestep, (u, vm) in solver.solve((0, tstop), dt):
        #plot(u, **plot_args)
        x = u.vector()[N/2-d]
        y = u.vector()[N/2+d]
        print u.vector().max(), x, y
        if x > threshold and (x-threshold)*(xp-threshold) < 0:
            # X cell has activated
            t0 = timestep[0]
        if y > threshold and (y-threshold)*(yp-threshold) < 0 and t0 != 0:
            # Y cell has activated
            t1 = timestep[0]
            # Calculate conduction velocity
            print "CV: ", 2*d*h/(t1-t0)*1e3
            break
        xp = x
        yp = y       


if __name__ == "__main__":
    BCL = 500
    dt = 0.1
    S2_time = 300
    tstop = 1500
    
    L = 250
    h = 0.5
    d = 10
    N = int(L/h)
    domain = RectangleMesh(0, 0, L, L, N, N)

    field_states = ["V"]
    field_parameters = ["stim_amplitude", "stim_offset"]
    save_fig = True

    plot_args = {'range_min': 0.,
                 'range_max': 1.,
                 'mode': 'color'}

    odes = ['hAM_KSMT_nSR', 'hAM_KSMT_cAF', 'FK_nSR', 'FK_cAF']
    thresholds = [0, 0, 0.5, 0.5]
    stim_amps = [-1410*9, -1410*9, -0.85, -0.85]
    D_list = [Constant(0.41), Constant(0.41), Constant(0.085), Constant(0.085)]

    S1_list = [Expression(("near(x[0],0)*amp","t"), amp=stim_amp, t=0) for stim_amp in stim_amps]
    S2_list = [Expression(("(x[1] < 0.5*L && x[0] < 0.5*L)*amp","t"),
                 L=L, t=S2_time+1, amp=stim_amp) for stim_amp in stim_amps]

    pa_list = [{'range_min':-80., 'mode':'color', 'range_max':40.}, {'range_min':-80., 'mode':'color', 'range_max':40.}, {'range_min':0., 'mode':'color', 'range_max':1.}, {'range_min':0., 'mode':'color', 'range_max':1.}]

    #plot(domain, interactive=True, mode='color')
    for i in range(2,len(odes),2):
        ode = odes[i]
        threshold = thresholds[i]
        D = as_tensor([[D_list[i], 0], [0, D_list[i]]])
        S1 = S1_list[i]
        S2 = S2_list[i]
        pa = pa_list[i]

        spiral(ode, domain, D, S1, S2, plot_args=pa, save_fig=True, threshold=threshold, BCL=500)
