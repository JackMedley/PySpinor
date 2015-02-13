#!/bin/python

# Import spinor library
from PySpinor import *

# Define momenta in the problem
pa = Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00, 7.500000e+02, incoming=True)  # Incoming gluon
pb = Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00,-7.500000e+02, incoming=True)  # Incoming gluon
p1 = Momenta(7.500000e+02, 1.663864e+02, 6.672462e+02,-2.993294e+02)		     # Outgoing ux
p2 = Momenta(7.500000e+02,-1.663864e+02,-6.672462e+02, 2.993294e+02)		     # Outgoing u

# Define physical spinors
u_a = Spinor(pa)
u_b = Spinor(pb)
u_1 = Spinor(p1)
u_2 = Spinor(p2)

pols = ('+', '+', '+', '+')

#
# t-channel
#

term_1 = Current(u_1.bar(), 'mu', u_a).dot(Current(u_2.bar(), 'mu', u_b))

A_t = - g_s ** 2 * term_1 / s(pa, p1)
print "A_t", A_t

#
# u-channel
#

# Sum over the unphysical spinor polarisations
term_2 = Current(u_2.bar(), 'mu', u_a).dot(Current(u_1.bar(), 'mu', u_b))

A_u = - g_s ** 2 * term_2 / s(pa, p2)
print "A_u", A_u

Sum = (A_t + A_u) * avgFactor() * 9.0

print "Sum for these polarisations = ",  Sum(pols)
