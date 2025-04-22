import numpy
import math

def monopole_power(r, c, rho, va, f):
    nfreq = f.shape[0]
    a = numpy.zeros([nfreq, 2])
    for i in range(nfreq):
        omega = 2 * math.pi * f[i]
        k = omega / c
        kr = k*r
        a[i,0] = f[i]
        a[i,1] = 2 * math.pi * r**2 * va**2 * rho * c * kr**2 / (1 + kr**2)
    numpy.savetxt("output.csv", a, delimiter=",")   

r = 0.5 # radius of sphere
f = numpy.linspace(10.0,1000.0,50)
c = 1500.0 # sound speed
rho = 1000.0 # density of fluid
va = 1.0

monopole_power(r, c, rho, va, f)
