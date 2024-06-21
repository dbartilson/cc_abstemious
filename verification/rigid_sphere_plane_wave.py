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

def rigid_sphere_plane_wave(r, c, f, loc, N = 60):
    nfp = loc.shape[0]
    nfreq = f.shape[0]
    if nfreq == 1:
        omega = 2 * math.pi * f[0]
        k = omega / c
        factors = numpy.full(N,0+0j)
        kr = k*r
        for n in range(N):
            factors[n] = 1j**n * (2*n+1) * jn_prime(n,kr) / hn1_prime(n,kr)
        p_scatt = numpy.full(nfp,0+0j)
        a = numpy.zeros([nfp, 4])
        for i in range(nfp):
            R = loc[i,0]
            theta = loc[i,1]
            z = numpy.cos(theta)
            kR = k*R
            for n in range(N):
                p_scatt[i] += factors[n] * legendre_fun(n, z) * hn1(n, kR)
                a[i,0] = R
                a[i,1] = theta
                a[i,2] = p_scatt[i].real
                a[i,3] = p_scatt[i].imag
        numpy.savetxt("output.csv", a, delimiter=",")   
    elif nfp == 1:
        R = loc[0,0]
        theta = loc[0,1]
        z = numpy.cos(theta)
        factors = numpy.full(N,0.0)
        p_scatt = numpy.full(nfreq,0+0j)
        a = numpy.zeros([nfreq, 3])
        for n in range(N):
            factors[n] = legendre_fun(n, z)
        for i in range(nfreq):
            omega = 2 * math.pi * f[i]
            k = omega / c
            kr = k*r
            kR = k*R
            for n in range(N):
                p_scatt[i] += factors[n] * 1j**n * (2*n+1) * jn_prime(n,kr) / hn1_prime(n,kr) * hn1(n, kR)
            a[i,0] = f[i]
            a[i,1] = p_scatt[i].real
            a[i,2] = p_scatt[i].imag
        numpy.savetxt("output.csv", a, delimiter=",")   

N = 80 # number of terms in sum
r = 0.5 # radius of sphere
R = 10.0 # radius of field points
nfp = 1 # number of field points
#f = numpy.logspace(math.log10(1.0),math.log10(1000.0),50) # frequency
f = numpy.linspace(10.0,1000.0,50)
c = 1500.0 # sound speed
# rho = 1000.0 # density of fluid

#theta = numpy.linspace(0, 2.0 * math.pi, num=nfp, endpoint=False)
theta = 0.0
loc = numpy.full([nfp,2],R)
loc[:,1] = theta

rigid_sphere_plane_wave(r, c, f, loc)
