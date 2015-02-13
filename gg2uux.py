#! /usr/bin/env python

# Import spinor library
from PySpinor import *

#
# Define momenta and spinors
#

pa = Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00, 7.500000e+02, incoming=True)  # Incoming gluon
pb = Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00,-7.500000e+02, incoming=True)  # Incoming gluon
p1 = Momenta(7.500000e+02, 1.663864e+02, 6.672462e+02,-2.993294e+02)		     # Outgoing ux
p2 = Momenta(7.500000e+02,-1.663864e+02,-6.672462e+02, 2.993294e+02)		     # Outgoing u

# Check the phase space point conserves momentum
assert Momenta.Conserved() is True, "ERROR! Momentum not conserved..."

u_a = Spinor(pa)
u_b = Spinor(pb)
u_1 = Spinor(p1)
v_1 = Spinor(p1, fermion=False)
u_2 = Spinor(p2)

# Select only one polarisation
pols = [('+', '-', '+', '-')]
pol  = pols[0]

# Define the gluon currents
Gluon_a = Gluon(pa, p1, 'nu')['+']
Gluon_b = Gluon(pb, p2, 'mu')['-']

#
# s-channel
#

term_1 = Current(u_2.bar(), 'mu', v_1).dot(pb - pa) * Gluon_b.dot(Gluon_a)
term_2 = Current(u_2.bar(), 'mu', v_1).dot(Gluon_a) * Gluon_b.dot(pa)
term_3 = Current(u_2.bar(), 'mu', v_1).dot(Gluon_b) * Gluon_a.dot(pb)

A_s = (term_1 + 2.0 * term_2 - 2.0 * term_3) / s(pa, pb)

print "   -> A_s = ", A_s[pol]

#
# t-channel
#

term_4  = Current(u_2.bar(), 'mu', u_a('+')).dot(Gluon_b) * Current(u_a.bar('+'), 'mu', v_1).dot(Gluon_a)
term_4 += Current(u_2.bar(), 'mu', u_a('-')).dot(Gluon_b) * Current(u_a.bar('-'), 'mu', v_1).dot(Gluon_a)

A_t = - i * term_4 / s(pa, p1)

print "   -> A_t = ", A_t[pol]

#
# u-channel
#

term_5  = Current(u_2.bar(), 'mu', u_a('+')).dot(Gluon_a) * Current(u_a.bar('+'), 'mu', v_1).dot(Gluon_b)
term_5 += Current(u_2.bar(), 'mu', u_a('-')).dot(Gluon_a) * Current(u_a.bar('-'), 'mu', v_1).dot(Gluon_b)
term_6  = Current(u_2.bar(), 'mu', u_2('+')).dot(Gluon_a) * Current(u_2.bar('+'), 'mu', v_1).dot(Gluon_b)
term_6 += Current(u_2.bar(), 'mu', u_2('-')).dot(Gluon_a) * Current(u_2.bar('-'), 'mu', v_1).dot(Gluon_b)

A_u = i * (term_5 - term_6) / s(pa, p2)

print "   -> A_u = ", A_u[pol]

#
# Form colour flow amplitudes for each polarisation
#

j0 = []
j1 = []

for pol in pols:

	j0.append((- i * A_s + A_u)[pol])
	j1.append((  i * A_s + A_t)[pol])

#
# Sum over the squared amplitudes with the colour factors and couplings
#

Sum = 0.0
for i, pol in enumerate(pols):

	print "\n   -> Matrix Element squared for polarisation", pol, " = ",

	temp = g_s ** 4 * (16.0 * abs2(j0[i]) + 16.0 * abs2(j1[i]) - 4.0 * (j0[i] * j1[i].conjugate()).real) / 3.0

	# Priunt indiviual amplitudes
	print temp

	# Add the contribution to the total sum
	Sum += temp

print "\n   -> Sum over all polarisation(s) = ",  Sum, "\n"


