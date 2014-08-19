from dolfin import *
from gotran import load_ode
from goss import DOLFINODESystemSolver, dolfin_jit
import numpy as np
import sys
sys.path.append("../0D/")
from steadystate import find_steadystate

def pacing_cv(ode, BCL_range, D, L, h, dt):
	"""
	Pace a 0D cell model for 50 cycles at given BCL, then
	simulate a wave in a 1D strand. 5 beats are initiated 
	at the left hand side of the strand, after the fifth 
	beat, the conduction velocity is measured.
	"""

	# Create domain
	N = int(L/h) # number of nodes - 1
	domain = IntervalMesh(N, 0, L)

	# Define functions
	V = FunctionSpace(domain, 'Lagrange', 1)
	u = TrialFunction(V)
	up = Function(V)
	v = TestFunction(V)

	# Assemble the mass matrix and the stiffness matrix
	M = assemble(u*v*dx)
	K = assemble(Constant(dt)*inner(D*grad(u), grad(v))*dx)
	A = M + K
	pde_solver = LUSolver(A)
	pde_solver.parameters["reuse_factorization"] = True

	# Set up ODESolver
	ode = load_ode(ode)
	compiled_ode = dolfin_jit(ode, field_states=["V"], field_parameters = [
							  "stim_period"])
	
	params = DOLFINODESystemSolver.default_parameters()
	params["scheme"] = "RL1"
	stim_field = interpolate(Expression("near(x[0],0)*sd", sd=1.), V)
	#compiled_ode.set_parameter("stim_duration", stim_field)
	BCL = 1000
	#compiled_ode.set_parameter("stim_period", BCL)

	ode_solver = DOLFINODESystemSolver(domain, compiled_ode, params=params)
	solver = ode_solver._ode_system_solvers[0]

	try:
		states = np.load("../data/steadystate_%s_BCL%d.npy" % (ode, BCL))
	except:
		states = find_steadystate(BCL, 50, dt, ode, plot_results=False)

	for i in range(N+1):
		solver.states(i)[:] = states
		up.vector()[i] = solver.states(i)[0] 

	#params = np.zeros((N+1)*2)
	solver.set_field_parameters(np.array([BCL]+[0]*N, dtype=np.float_))
	#stim_field = interpolate(Expression("0"), V)
	#compiled_ode.set_parameter("stim_duration", stim_field)



	t = 0.
	tstop = 1e6
	while t < tstop:
		# Step ODE solver
		ode_solver.step((t, t+dt), up)

		# Assemble RHS and solve pde
		b = M*up.vector()
		pde_solver.solve(up.vector(), b)

		# Plot solution
		plot(up, range_min=-80.0, range_max=40.0)


		t += dt

ode = "hAM_KSMT_nSR"
BCL_range=[1000]
D = 5e-7
L = 0.1
h = 0.0005
dt = 0.05

pacing_cv(ode, BCL_range, D, L, h, dt)


'''
	for BCL in BCL_range():
		try:
			states = np.load("../data/steadystate_%s_BCL%d.npy" % (ode, BCL))
		except:
			states = find_steadystate(BCL, 50, dt, ode, plot_results=False)

		compiled_ode.set_parameter("stim_period", BCL)
		

		BCL = 1000
    num_of_beats = 50
    dt = 0.01
    ode = "hAM_KSMT_nSR"

    find_steadystate(BCL, num_of_beats, dt, ode, plot_results=True)




# Parameters
D = Constant(1e-4) # Diffusion constant
stim_duration = 1
stim_amp = -1410*2
BCL = 400

# Time resolution
dt = 0.01
tstop = 6*BCL


# ODESolver
ode = load_ode("hAM_KSMT_nSR.ode")
ode_compiled = dolfin_jit(ode, field_states=["V"], field_parameters=["stim_duration", "stim_period", "stim_amplitude"])
params = DOLFINODESystemSolver.default_parameters()
# Allowed schemes: BasicImplicitEuler, ESDIRK23a, ESDIRK4032, ExplicitEuler, GRL1, GRL2, ImplicitEuler, RK2, RK4, RKF32, RL1, RL2, ThetaSolver
params["scheme"] = "RL1" 
stim_field = interpolate(Expression("near(x[0],0)*sd", sd=stim_duration), V)
ode_compiled.set_parameter("stim_duration", stim_field)
ode_compiled.set_parameter("stim_period", BCL)
ode_compiled.set_parameter("stim_amplitude", stim_amp)

ode_solver = DOLFINODESystemSolver(domain, ode_compiled, params=params)
#solver.reset_default()

solver = ode_solver._ode_system_solvers[0]

states  = np.load("steadystate.npy")
for i in range(N+1):
	solver.states(i)[:] = states
	up.vector()[i] = solver.states(i)[0] 
	print up.vector()[i]

#ode_solver.reset_default()
#ode_solver.states(1)

plot(up, range_min=-80.0, range_max=40.0)

# Iterate 
t = 0.
xp = yp =-1;
t1 = t2 = 0;
while float(t) < tstop:
	dt = 0.01 if (float(t) % BCL < 100) else 0.1

	# Step ODE solver
	ode_solver.step((t, t+dt), up)
	
	# Assemble rhs
	b = M*up.vector()

	# Solve
	pde_solver.solve(up.vector(), b)

	x = up.vector().array()[N/2]
	y = up.vector().array()[(N/2)+1]

	if x*xp < 0 and xp < 0:
		t1 = t
	if y*yp < 0 and yp < 0:
		t2 = t
	if t1 and t2:
		print "C_v = ", 10**6*float(h)/(t2-t1)
		t1 = t2 = 0
	xp = x
	yp = y

 	# Plot solution
	plot(up, range_min=-80.0, range_max=40.0)

	t += dt
'''