from gotran.codegeneration import compile_module
from gotran.common.options import parameters
from gotran.model.loadmodel import load_ode
import numpy as np

class CellModel:
    """
    Class for compiling a gotran ODE cell model and 
    perform pacing experiments on it.
    """

    def __init__(self, ode, solver="rush_larsen",
                            monitored=[],
                            odepath="../ode/", 
                            scpath="../data/steadycycles"):
        """
        Compiles a gotran ode and stores the module and forward method.
        """
        generation = parameters.generation.copy()

        # We don't want to compile unessecery features
        for what_not in ["componentwise_rhs_evaluation",
                         "forward_backward_subst",
                         "linearized_rhs_evaluation",
                         "lu_factorization",
                         "jacobian"]:
            generation.functions[what_not].generate = False

        # Select what numerical scheme we want the solver to use
        generation.solvers[solver].generate = True

        if isinstance(ode, str):
            ode = load_ode(odepath+ode)
        
        self.module = compile_module(ode, "C", monitored, generation)
        self.forward = getattr(self.module, "forward_"+solver)
        self.solver = solver
        self.odepath = odepath
        self.scpath = scpath
        self.index = {}

    def get_state_indices(*states):
        for state in states:
            self.index[state] = self.module.state_indices[state]

    def get_parameter_indices(*params):
        for param in params:
            self.index[param] = self.module.parameter_indices[param]

    def set_dt(self, dt):
        self.dt = dt

    def set_odepath(self, odepath):
        self.odepath = odepath

    def set_scpath(self, scpath):
        self.scpath = scpath

cell = CellModel('FK_nSR')
cell.get_state_indices('V')
cell.get_paramer_indices('stim_period')


#     def pace()


# def pacing_S1S2(ode, S1, S2_range, dt, threshold=0.1, plot=False,
#                 odepath='../ode/', scpath="../data/steadycycles/"):
#     """
#     Calculate a restitution curve using S1-S2 pacing.
#     The ODE model is paced for 50 beats at every BCL value,
#     and the APD and DI are measured for the final beat.
#     """

#     # Load in inital state after S1 pacing
#     try:
#         init_states = np.load(scpath+"%s_BCL%d.npy" % (ode, S1))
#     except:
#         print "Steady cycle at S1=%d for ODE model: %s not found." % (S1, ode)
#         print "Pacing 0D cell model to find it, this may take a minute."
#         states = find_steadycycle([module, forward, ode], BCL, dt, odepath=odepath, scpath=scpath)
#         print "Steady cycle found, proceeding to do S1-S2 pacing."

#     # Get state and parameter indices
#     model_params = module.init_parameter_values()
#     index = {}
#     index['V'] = module.state_indices('V')
#     index['BCL'] = module.parameter_indices('stim_period')
    
#     def pulse(S2):
#         """Pulse the steady cycle with a single S2 pulse and measure APD"""
#         states = init_states.copy()
#         model_params[index['BCL']] = S2

#         while t <= S2:
#             forward(states, t, dt, model_params)
#             t += dt

#         V = [states[index['V']]]
#         while t <= 2*S2:
#             forward(states, t, dt, model_params)
#             t += dt
#             V.append(states[index['V']])

#         # Extract APD
#         V = np.array(V)
#         APD = len(V[V>threshold])*dt
#         DI = S2 - APD
#         return APD, DI

#     for S2 in S2_range:
#         APD, DI = pulse(S2)
#         print "S2: %g,\t APD: %g,\t DI: %g" % (S2, APD, DI)

# if __name__ == '__main__':
#     ### Example of use

#     ode = 'hAM_KSMT_nSR'
#     ode = 'hAM_KSMT_cAF'
#     ode = 'FK_nSR'
#     ode = 'FK_cAF'

#     BCL_range = range(600, 495, -5)
#     dt = 0.01
#     S1 = 1000
#     S2_range = range(500, 295, -5)

#     pacing_S1S2(ode, S1, S2_range, dt)
