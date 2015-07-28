from TissueSquare import *

def random_dead(D, p):
    return Expression((("(rand() > p)*D", "0"), ("0", "(rand() > p)*D")), D=D, p=p)

def fib_heart(scale):
    fib_scale = Constant(scale)
    absx = 'std::abs(4*(2*x[0]/L-(1-c))/(2*c)-2)'
    shape = "(1 - (x[1] < (1-c)/2.0*L + c*L/4.0*(3+sqrt(1-(%s-1)*(%s-1))) && x[1] > (1-c)/2.0*L + c*L/4.0*(3-3*sqrt(1-sqrt(%s/2)))))" % (absx, absx, absx)
    D_scale = Constant(0.38)
    fib_shape = Expression((("D*"+shape, "0"), ("0", "D*"+shape)), L=L, c=fib_scale, D=D_scale)
    return fib_shape

ode = 'FK_nSR'
# D = random_dead(0.26, 0.65)
amp = -0.8
plot = {'range_min':0.0, 'range_max':1.0, 'mode':'color', 'interactive':False}
BCL = 500
L = 50
D = random_dead(0.26, 0.8)


solver = TissueSquare(ode, D, L=L, stim_amp=amp, plot_args=plot)
solver.spiral_wave(BCL, S2time=200, tstop=1000, dt=0.1, savefig=True)



'''
uplot = u.compute_vertex_values()
uplot = uplot.reshape(sqrt(uplot.size), sqrt(uplot.size))
plt.imshow(uplot, interpolation='none')
plt.show()
'''