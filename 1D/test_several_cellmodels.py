from dolfin import *
from gotran import load_ode
from goss import dolfin_jit
from cbcbeat.cardiacmodels import CardiacModel
from cbcbeat.gossplittingsolver import GOSSplittingSolver
from cbcbeat.markerwisefield import Markerwise
import numpy as np
import matplotlib.pyplot as plt

set_log_level(ERROR)
# Load states into 1D model

def GOSSparams():
    params = GOSSplittingSolver.default_parameters()
    params["pde_solver"] = "monodomain"
    params["MonodomainSolver"]["linear_solver_type"] = "iterative"
    params["MonodomainSolver"]["theta"] = 1.0
    params["ode_solver"]["solver"] = "RL1"
    params["apply_stimulus_current_to_pde"] = False
    return params

domain = IntervalMesh(50, 0, 10)

# Load the cell model from file
ode1 = "hAM_KSMT_nSR"
ode2 = "hAM_KSMT_cAF"

lode1 = load_ode(ode1)
lode2 = load_ode(ode2)

cellmodel1 = dolfin_jit(lode1, field_states=['V'],
                       field_parameters=['stim_period', 'stim_amplitude'])

cellmodel2 = dolfin_jit(lode2, field_states=['V'],
                       field_parameters=['stim_period', 'stim_amplitude'])

markers = MeshFunction("size_t", domain, 0, 1)

passive = CompiledSubDomain("x[0] > 5 && 6 > x[0]")
passive.mark(markers, 2)

cellmodels = Markerwise((cellmodel1, cellmodel2), (1,2), markers)

M_i = Constant(0.34) #as_tensor([[Constant(0.34), 0],[0, Constant(0.34)]])
M_e = Constant(0.34) #as_tensor([[Constant(0.34), 0],[0, Constant(0.34)]])
t = Constant(0.0)

# Create the CardiacModel for the given mesh and cell model
tissue = CardiacModel(domain, t, M_i, M_e, cellmodels)

# Create the solver and extract the subsolvers
solver = GOSSplittingSolver(tissue, GOSSparams())
dolfin_solver = solver.ode_solver # Solves spatial PDE
ode_solver1 = dolfin_solver._ode_system_solvers[1] # Solves ODE system
ode_solver2 = dolfin_solver._ode_system_solvers[2] # Solves ODE system

# Initialize cell by pacing 0D cell model
states1 = np.load("../data/steadycycles/" + "%s_BCL%d.npy" % (ode1, 300))
states2 = np.load("../data/steadycycles/" + "%s_BCL%d.npy" % (ode2, 300))

# Set the stimulus parameter field
N=50
L=10
h=L/50.

stim_field1 = np.zeros(2*ode_solver1.num_nodes(), dtype=np.float_)
stim_field2 = np.zeros(2*ode_solver2.num_nodes(), dtype=np.float_)

stim_field1[::2] = 300
stim_field2[::2] = 300

stim_field1[-1] = -1410*12
#stim_field[] = -1410*8

ode_solver1.set_field_parameters(stim_field1)
ode_solver2.set_field_parameters(stim_field2)

# Load states into 1D model
(u, um) = solver.solution_fields()

for node in range(ode_solver1.num_nodes()):
    ode_solver1.states(node)[:] = states1.copy()

for node in range(ode_solver2.num_nodes()):
    ode_solver2.states(node)[:] = states2.copy()

u.vector()[:] = states1[0]
 
for timestep, (u, vm) in solver.solve((0, 1000), 0.1):
    plot(u, range_min=-80.0, range_max=40.0)