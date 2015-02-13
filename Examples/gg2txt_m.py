#! /usr/bin/env python

# Import spinor library
from PySpinor import *

#
# Define momenta in the problem
#

pa = Momenta(7.500000e+02,  0.000000e+00,  0.000000e+00,  7.500000e+02, incoming=True)  # Incoming gluon
pb = Momenta(7.500000e+02,  0.000000e+00,  0.000000e+00, -7.500000e+02, incoming=True)  # Incoming gluon
p1 = Momenta(7.500000e+02,  1.618995e+02,  6.492524e+02, -2.912573e+02)		             # Outgoing tx
p2 = Momenta(7.500000e+02, -1.618995e+02, -6.492524e+02,  2.912573e+02)		             # Outgoing t

# Check the phase space point conserves momentum
assert Momenta.Conserved() is True, "ERROR! Momentum not conserved..."

m = p1.mass()

# Now things are more complicated since the tops have a non-negligible mass
# Start by decomposing the top momenta into two basis momenta (each)
k1, k2 = p1.decompose()
k3, k4 = p2.decompose()

# Spinors
u_a = Spinor(pa)
u_b = Spinor(pb)
u_1 = Spinor(p1)
u_2 = Spinor(p2)
v_1 = Spinor(p1, fermion=False)
v_2 = Spinor(p2, fermion=False)

# Dirac Equation Checks
# print p1.slashed() * u_1('+')
# print m * u_1('+'), "\n"
# print p1.slashed() * v_1('+')
# print m * v_1('+'), "\n"
# print p2.slashed() * u_2('+')
# print m * u_2('+'), "\n"
# print p2.slashed() * v_2('+')
# print m * v_2('+'), "\n"

# Spinor completeness checks
# print u_1('+') * u_1.bar('+') + u_1('-') * u_1.bar('-')
# print "\n"
# print p1.slashed() + m * I4
# print "\n"
# print v_1('+') * v_1.bar('+') + v_1('-') * v_1.bar('-')
# print "\n"
# print p1.slashed() - m * I4
# print "\n"

# Pick polarisations
polarisations = [('+', '-', '+', '-')]
pol = polarisations[0]

# Choose two massless reference momenta:
# Changes in this should be equivalent to a gauge transformation
ra = k2
rb = k3

# Define the gluon currents using these reference momenta
Gluon_a = Gluon(pa, ra, 'nu')['+']
Gluon_b = Gluon(pb, rb, 'mu')['-']

#
# s
#

term_1 = Current(u_2.bar('-'), 'mu', v_1('+')).dot(pb - pa) * Gluon_b.dot(Gluon_a)
term_2 = Current(u_2.bar('-'), 'mu', v_1('+')).dot(Gluon_a) * Gluon_b.dot(pa)
term_3 = Current(u_2.bar('-'), 'mu', v_1('+')).dot(Gluon_b) * Gluon_a.dot(pb)

A_s = (term_1 + 2.0 * term_2 - 2.0 * term_3) / s(pa, pb)

print "   -> A_s = ", A_s[pol]

#
# t
#

# With the antispinor completion relation
term_1 = Current(u_2.bar('-'), 'mu', u_a('+')).dot(Gluon_b) * \
         Current(u_a.bar('+'), 'nu', v_1('+')).dot(Gluon_a) + \
         Current(u_2.bar('-'), 'mu', u_a('-')).dot(Gluon_b) * \
         Current(u_a.bar('-'), 'nu', v_1('+')).dot(Gluon_a)

term_2 = Current(u_2.bar('-'), 'mu', v_1('+')).dot(Gluon_b) * \
         Current(v_1.bar('+'), 'nu', v_1('+')).dot(Gluon_a) + \
         Current(u_2.bar('-'), 'mu', v_1('-')).dot(Gluon_b) * \
         Current(v_1.bar('-'), 'nu', v_1('+')).dot(Gluon_a)

A_t = - i * (term_1 - term_2) / s(pa, p1)

print "   -> A_t = ", A_t[pol]

#
# u
#

term_1 = Current(u_2.bar('-'), 'nu', u_a('+')).dot(Gluon_a) * \
         Current(u_a.bar('+'), 'mu', v_1('+')).dot(Gluon_b) + \
         Current(u_2.bar('-'), 'nu', u_a('-')).dot(Gluon_a) * \
         Current(u_a.bar('-'), 'mu', v_1('+')).dot(Gluon_b)

term_2 = Current(u_2.bar('-'), 'nu', u_2('+')).dot(Gluon_a) * \
         Current(u_2.bar('+'), 'mu', v_1('+')).dot(Gluon_b) + \
         Current(u_2.bar('-'), 'nu', u_2('-')).dot(Gluon_a) * \
         Current(u_2.bar('-'), 'mu', v_1('+')).dot(Gluon_b)

A_u = i * (term_1 - term_2) / s(pa, p2)

print "   -> A_u = ", A_u[pol]

# 
# Form colour flow amplitudes for each polarisation
#

j0 = []
j1 = []

for pol in polarisations:

   j0.append((- i * A_s + A_u)[pol])
   j1.append((  i * A_s + A_t)[pol])

#
# Sum over the squared amplitudes with the colour factors and couplings
#

Sum = 0.0
itr = 0
for i, pol in enumerate(polarisations):

	print "\n   -> Matrix Element squared for polarisation", polarisations[i], " = ",

	temp = g_s ** 4 * (16.0 * abs2(j0[i]) + 16.0 * abs2(j1[i]) - 4.0 * (j0[i] * j1[i].conjugate()).real) / 3.0

	# Priunt indiviual amplitudes
	print temp

print "\n"
if itr > 1:

	print "   -> Sum over all polarisation(s) = ",  Sum, "\n"
