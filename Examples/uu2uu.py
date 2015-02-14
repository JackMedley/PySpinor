#!/bin/python

# Import spinor library
import PySpinor as pys
from   PySpinor import Common
from   PySpinor.Common   import *

# Define momenta in the problem
pa = pys.Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00, 7.500000e+02, label='a', incoming=True)  # Incoming gluon
pb = pys.Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00,-7.500000e+02, label='b', incoming=True)  # Incoming gluon
p1 = pys.Momenta(7.500000e+02, 1.663864e+02, 6.672462e+02,-2.993294e+02, label='1')		    # Outgoing ux
p2 = pys.Momenta(7.500000e+02,-1.663864e+02,-6.672462e+02, 2.993294e+02, label='2')		    # Outgoing u

# Define physical spinors
u_a = pys.Spinor(pa)
u_b = pys.Spinor(pb)
u_1 = pys.Spinor(p1)
u_2 = pys.Spinor(p2)

print u_a.momentum.__dict__

pols = ('+', '+', '+', '+')

#
# t-channel
#

term_1 = pys.Current(u_1.bar(), 'mu', u_a).dot(pys.Current(u_2.bar(), 'mu', u_b))

A_t = - g_s ** 2 * term_1 / pys.s(pa, p1)

#
# u-channel
#

# Sum over the unphysical spinor polarisations
term_2 = pys.Current(u_2.bar(), 'mu', u_a).dot(pys.Current(u_1.bar(), 'mu', u_b))

A_u = - g_s ** 2 * term_2 / pys.s(pa, p2)

Sum = (A_t + A_u) * pys.avgFactor() * 9.0

print "Sum for these polarisations = ",  Sum(pols)

# Test of the (non-automated) printing of latex
printOut = pys.spinorString(u_1.bar(), 'mu', u_a, u_2.bar(), 'mu', u_b, prefactor='-\\frac{g_s^2}{s_{a1}}')

printOut.addString()

printOut.getLaTeX('Test')

printOut.compileLaTeX()