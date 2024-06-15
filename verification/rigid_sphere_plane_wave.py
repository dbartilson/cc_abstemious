import scipy
import numpy
import math
import cmath

# derivative of spherical besselj
def jn_prime(n, z):
    return (-scipy.special.spherical_jn(n+1,z) + (n/z)*scipy.special.spherical_jn(n,z))

# spherical hankel function of first kind
def hn1(n, z):
    return scipy.special.spherical_jn(n,z) + 1j*scipy.special.spherical_yn(n,z)

# derivative of spherical hankel function of first kind
def hn1_prime(n, z):
    return (-hn1(n+1,z) + (n/z)*hn1(n,z))

# legendre function (i.e., legendre polynomal evaluated at ordinate)
def legendre_fun(n, z):
    p = scipy.special.legendre(n)
    return p(z)

N = 80 # number of terms in sum
r = 0.5 # radius of sphere
R = 10.0 # radius of field points
nfp = 100 # number of field points
f = 10.0 # frequency
c = 1500.0 # sound speed
rho = 1000.0 # density of fluid

omega = 2 * math.pi * f
k = omega / c

theta = numpy.linspace(0, 2.0 * math.pi, num=nfp, endpoint=False)

factors = numpy.full(N,0+0j)
kr = k*r
for n in range(N):
    factors[n] = 1j**n * (2*n+1) * jn_prime(n,kr) / hn1_prime(n,kr) 

p_scatt = numpy.full(nfp,0+0j)
a = numpy.zeros([nfp, 3])
kR = k*R
for i in range(nfp):
    z = numpy.cos(theta[i])
    for n in range(N):
        p_scatt[i] += factors[n] * legendre_fun(n, z) * hn1(n, kR)
        a[i,0] = theta[i]
        a[i,1] = p_scatt[i].real
        a[i,2] = p_scatt[i].imag

numpy.savetxt("output.csv", a, delimiter=",")