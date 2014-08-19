from dolfin import *
from gotran import load_ode
from goss import dolfin_jit
from beatadjoint.cardiacmodels import CardiacModel
from beatadjoint.gossplittingsolver import GOSSplittingSolver
import numpy as np

# Turn of dolfin printing
set_log_level(ERROR)

class Monodomain2D:
    def __init__(self, ode):
        # Create the domain
        domain = RectangleMesh(0, 0, L, L, N, N)

        # Load the cell model from file
        cellmodel = load_ode(ode)
        cellmodel = dolfin_jit(cellmodel, field_states, field_parameters)

        # Create the CardiacModel for the given domain and cell model
        heart = CardiacModel(domain, Constant(0.0), D, None, cellmodel)

        # Create the solver
        solver = GOSSplittingSolver(heart, self.GOSSparams())
        
        # Get the solution fields and subsolvers
        dolfin_solver = solver.ode_solver
        ode_solver = dolfin_solver._ode_system_solvers[0]

        self.solver = solver
        self.ode = ode
        self.dolfin_solver = dolfin_solver
        self.ode_solver = ode_solver
        self.domain = domain

    def GOSSparams(self):
        params = GOSSplittingSolver.default_parameters()
        params["pde_solver"] = "monodomain"
        params["MonodomainSolver"]["linear_solver_type"] = "iterative"
        params["MonodomainSolver"]["theta"] = 1.0
        params["ode_solver"]["scheme"] = "RL1"
        params["apply_stimulus_current_to_pde"] = False
        return params

    def diffusion_constant(self):
        Dx = Constant(0.4)
        Dy = Constant(0.4)
        D = as_tensor([[Dx, 0], [0, Dy]])
        return D

    def domain():
        L = 15   # size of mesh, [mm]
        h = 0.1  # spatial resolutiofn, [mm]
        dt = 0.1 # temporal resolution, [ms]

        # Create dolfin mesh
        N = int(L/h) # number of nodes in each direction
        domain = RectangleMesh(0, 0, L, L, N, N)
        return domain, L

    def planar_wave(self):
        '''
        Finds the conduction velocity of a planar wave. First paces the 0D ode
        model for 50 cycles at a BCL of 1000, then initiates a single planar
        wave, measuring the conduction velocity half-way in the domain.
        '''
        # Import function for pacing 0D cell model to find quasi steady-state
        import sys
        sys.path.append("../0D/")
        from steadystate import find_steadystate
        
        solver, ode_solver = self.solver, self.ode_solver

        # Pace 0D cell model and find quasi-steady state
        try:
            steadystate = np.load("../data/steadystate_%s_BCL%d.npy"%(ode,1000))
        except:
            BCL = 1000
            print "Pacing 0D model for BCL %g..." % BCL,
            steadystate = find_steadystate(BCL, 50, dt, ode, plot_results=False)
            print "done."

        # Get the solution fields
        (u, um) = solver.solution_fields()

        # Load steadystate into 2D model
        for i in range((N+1)*(N+1)):
            ode_solver.states(i)[:] = steadystate
            u.vector()[i] = steadystate[0]
        
        # Define the planar wave stimulus
        S = Expression(("near(x[0],0)*amp", "t"), amp=stim_amp, t=0)
        V = VectorFunctionSpace(self.domain, 'CG', 1)
        S = interpolate(S, V).vector().array()

        # Apply the stimulus
        ode_solver.set_field_parameters(S)

        xp = 0; yp = 0;
        for timestep, (u, vm) in solver.solve((0, tstop), dt):
            plot(u, **plot_args)
            print u.vector().max()
            x = u.vector()[N/2-d]
            y = u.vector()[N/2+d]
            if x > 0 and (x*xp) < 0:
                # X cell has activated
                t0 = timestep[0]
            if y > 0 and (y*yp) < 0 and t0 != 0:
                # Y cell has activated
                t1 = timestep[0]
                # Calculate conduction velocity
                print "CV: ", 2*d*h/(t1-t0)*1e3
                #break
            xp = x
            yp = y

    def spiral_wave(self):
        solver, ode_solver = self.solver, self.ode_solver
    
        # Define stimulus fields
        S1 = Expression(("near(x[0],0)*amp","t"), amp=stim_amp, t=0)
        S2 = Expression(("(x[1] < 0.5*L && x[0] < 0.5*L)*amp","t"), L=L, t=S2_time+1, amp=stim_amp)
        nostim = Expression(("0","0"))

        # Calculate parameter fields 
        V = VectorFunctionSpace(self.domain, 'CG', 1)
        S1 = interpolate(S1, V).vector().array()
        S2 = interpolate(S2, V).vector().array()
        nostim = interpolate(nostim, V).vector().array()
        
        # Pace 0D cell model and find quasi-steady state
        try:
            steadystate = np.load("../data/steadystate_%s_BCL%d.npy"%(ode, BCL))
        except:
            print "Pacing 0D model for BCL %g..." % BCL,
            steadystate = find_steadystate(BCL, 50, dt, ode, plot_results=False)
            print "done."

         # Get the solution fields
        (u, um) = solver.solution_fields()

        # Load steadystate into 2D model
        for i in range((N+1)*(N+1)):
            ode_solver.states(i)[:] = steadystate
            u.vector()[i] = steadystate[0]

        cnt = 0
        plot(u, interactive=True, **plot_args)

        # Stimulate S1 pulse
        ode_solver.set_field_parameters(S1)
        for timestep, (u, vm) in solver.solve((0, S2_time), dt):
            fig = plot(u, interactive=False, **plot_args)

            if save_fig and cnt % 5 == 0:
                padded_index = '%08d' % cnt
                fig.write_png('../tmp/fib_%s' % self.ode + padded_index)
            cnt += 1

        ode_solver.set_field_parameters(S2)

        # Stimulate S2 pulse
        for timestep, (u, vm) in solver.solve((S2_time, S2_time+100), dt):
            fig = plot(u, interactive=False, **plot_args)

            if save_fig and cnt % 5 == 0:    
                padded_index = '%08d' % cnt
                fig.write_png('../tmp/fib_%s' % self.ode + padded_index)
            cnt += 1

        # Turn of all stimulus and iterate until simulation stop
        ode_solver.set_field_parameters(nostim)
        for timestep, (u, vm) in solver.solve((S2_time+100, tstop), dt):
            fig = plot(u, interactive=False, **plot_args)

            if save_fig and cnt % 5 == 0:
                padded_index = '%08d' % cnt
                fig.write_png('../tmp/fib_%s' % self.ode + padded_index)
            cnt += 1


if __name__ == "__main__":
    # Parameters
    BCL = 500

    dt = 0.1
    S2_time = 300
    tstop = 1500
    L = 250
    h = 0.5
    d = 10
    N = int(L/h)
    stim_amp = -1410*9
    field_states = ["V"]
    field_parameters = ["stim_amplitude", "stim_offset"]
    axis = [-80., 40.]
    save_fig = True

    plot_args = {'range_min': -80.,
                 'range_max': 40.,
                 'mode': 'color',
                 #'window_width': 600,
                 #'window_height': 400
                 }

    fib_scale = Constant(0.20)
    absx = 'std::abs(4*(2*x[0]/L-(1-c))/(2*c)-2)'
    shape = "(1 - (x[1] < (1-c)/2.0*L + c*L/4.0*(3+sqrt(1-(%s-1)*(%s-1))) && x[1] > (1-c)/2.0*L + c*L/4.0*(3-3*sqrt(1-sqrt(%s/2)))))" % (absx, absx, absx)
    D_scale = Constant(0.38)
    fib_shape = Expression((("D*"+shape, "0"), ("0", "D*"+shape)), L=L, c=fib_scale, D=D_scale)
    D = fib_shape

    # Produce spiral wave
    odes = ['hAM_KSMT_cAF']
    for ode in odes:
        solver = Monodomain2D(ode)
        solver.spiral_wave()

