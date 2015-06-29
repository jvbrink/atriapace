from pylab import *

Vo = 0.0
Vm = 1.0
Vna = 0.23
Vc = 0.33215242
Vv = 0.67095874
Vw = 0.27314403
Vd = 0.33471235
tvm = 122.36631378
tvp = 2.00966394
twm  = 210.11018188
twp = 627.77222072
tsp = 0.38651189
tsm = 0.27632535
Vcsi = 0.23248197
xk  = 5.81406155
td = 0.04237912
to  = 22.18688450
tsoa = 47.70819517
tsob = 2.44640465
Vso = 0.60666388
xtso = 5.85023408
tsi = 46.51758570
D = 0.00091753
tvmm = 1300.32077428

V = linspace(0,1,1001)

v = 0

def Jfi(v):
    return -v/td * (Vc < V)*(1-V)*(V-Vc)

plot(V, Jfi(0))
plot(V, Jfi(0.25))
plot(V, Jfi(0.5))
plot(V, Jfi(0.75))
plot(V, Jfi(1))
show()