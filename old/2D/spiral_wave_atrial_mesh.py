from dolfin import *
from gotran import load_ode
from goss import dolfin_jit
from beatadjoint.cardiacmodels import CardiacModel
from beatadjoint.gossplittingsolver import GOSSplittingSolver
import numpy as np

# Turn of dolfin printing
set_log_level(ERROR)

BCL = 200
Dx = Constant(0.005) #38)
Dy = Constant(0.005) #38)
Dz = Constant(0.005) #38)
D = as_tensor([[Dx, 0, 0], [0, Dy, 0], [0, 0, Dz]])
dt = 0.1
S2_time = 130
stim_period = 80
tstop = 1000
stim_amp = -0.3
field_states = ["V"]
field_parameters = ["stim_amplitude", "stim_offset", "stim_period"]
axis = [-80., 40.]
save_fig = True

plot_args = {'range_min': 0.,
             'range_max': 1.,
             'mode': 'color',
             #'window_width': 1366,
             #'window_height': 744
             }

ode = 'FK_cAF'
domain = Mesh('atrial_mesh.xml')

# Load the cell model from file
cellmodel = load_ode(ode)
cellmodel = dolfin_jit(cellmodel, field_states, field_parameters)

def GOSSparams():
    params = GOSSplittingSolver.default_parameters()
    params["pde_solver"] = "monodomain"
    params["MonodomainSolver"]["linear_solver_type"] = "iterative"
    params["MonodomainSolver"]["theta"] = 1.0
    params["ode_solver"]["scheme"] = "RL1"
    params["apply_stimulus_current_to_pde"] = False
    return params

# Create the CardiacModel for the given domain and cell model
heart = CardiacModel(domain, Constant(0.0), D, None, cellmodel)

# Create the solver
solver = GOSSplittingSolver(heart, GOSSparams())
    
# Get the solution fields and subsolvers
dolfin_solver = solver.ode_solver
ode_solver = dolfin_solver._ode_system_solvers[0]

# Define stimulus fields
S1 = ("(x[2] > 15.1)*amp", "t", "10000")
#S2 = ("(x[1] > 4.75)*amp", "t", "10000")
S2 = ("(x[0] > 0.0 && x[0] < 2.0 && x[2] > 10.0)*amp", "t", "stim_period")
#("(x[0] > 4.8 && x[0] < 5.3 && x[1] > 4.8 && x[1] < 5.3)*amp", "t", "10000")
#("(x[2] > 14 && x[1] < 4)*amp", "t", "10000")
S1 = Expression(S1, amp=stim_amp, t=0)
S2 = Expression(S2, amp=stim_amp, t=S2_time+1, stim_period=stim_period)
nostim = Expression(("0","0","10000"))

# Calculate parameter fields 
V = VectorFunctionSpace(domain, 'CG', 1)
S1 = interpolate(S1, V).vector().array()
S2 = interpolate(S2, V).vector().array()
nostim = interpolate(nostim, V).vector().array()
    
# Pace 0D cell model and find quasi-steady state
try:
    steadystate = np.load("../data/steadystates/steadystate_%s_BCL%d.npy"%(ode, BCL))
except:
    print "Pacing 0D model for BCL %g..." % BCL,
    steadystate = find_steadystate(BCL, 50, dt, ode, plot_results=False)
    print "done."

# Get the solution fields
(u, um) = solver.solution_fields()

# Load steadystate into 2D model
for i in range(domain.coordinates().shape[0]):
    ode_solver.states(i)[:] = steadystate
    u.vector()[i] = steadystate[0]

cnt = 0
# Stimulate S1 pulse
plot(u, interactive=True, **plot_args)

ode_solver.set_field_parameters(S1)
for timestep, (u, vm) in solver.solve((0, 10), dt):
    fig = plot(u, interactive=False, **plot_args)

    if save_fig and cnt % 10 == 0:
        padded_index = '%08d' % cnt
        fig.write_png('../tmp/atrial_spiral_%s' % ode + padded_index)
    cnt += 1

ode_solver.set_field_parameters(S2)

# Stimulate S2 pulse
for timestep, (u, vm) in solver.solve((10, S2_time+5*stim_period), dt):
    fig = plot(u, interactive=False, **plot_args)

    if save_fig and cnt % 10 == 0:
        padded_index = '%08d' % cnt
        fig.write_png('../tmp/atrial_spiral_%s' % ode + padded_index)
    cnt += 1

# Turn of all stimulus and iterate until simulation stop
ode_solver.set_field_parameters(nostim)
for timestep, (u, vm) in solver.solve((S2_time+5*stim_period, tstop), dt):
    fig = plot(u, interactive=False, **plot_args)

    if save_fig and cnt % 10 == 0:
        padded_index = '%08d' % cnt
        fig.write_png('../tmp/atrial_spiral_%s' % ode + padded_index)
    cnt += 1

