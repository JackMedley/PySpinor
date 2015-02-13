#! /usr/bin/env python

# Import spinor library
import PySpinor as pys

#
# Define momenta and spinors
#

pa = pys.Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00, 7.500000e+02, incoming=True)  # Incoming gluon
pb = pys.Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00,-7.500000e+02, incoming=True)  # Incoming gluon
p1 = pys.Momenta(7.500000e+02, 1.663864e+02, 6.672462e+02,-2.993294e+02)		     # Outgoing ux
p2 = pys.Momenta(7.500000e+02,-1.663864e+02,-6.672462e+02, 2.993294e+02)		     # Outgoing u

# Check the phase space point conserves momentum
assert pys.Momenta.Conserved() is True, "ERROR! Momentum not conserved..."

u_a = pys.Spinor(pa)
u_b = pys.Spinor(pb)
u_1 = pys.Spinor(p1)
v_1 = pys.Spinor(p1, fermion=False)
u_2 = pys.Spinor(p2)

# Select only one polarisation
pols = [('+', '-', '+', '-')]
pol  = pols[0]

# Define the gluon currents
Gluon_a = pys.Gluon(pa, p1, 'nu')['+']
Gluon_b = pys.Gluon(pb, p2, 'mu')['-']

#
# s-channel
#

term_1 = pys.Current(u_2.bar(), 'mu', v_1).dot(pb - pa) * Gluon_b.dot(Gluon_a)
term_2 = pys.Current(u_2.bar(), 'mu', v_1).dot(Gluon_a) * Gluon_b.dot(pa)
term_3 = pys.Current(u_2.bar(), 'mu', v_1).dot(Gluon_b) * Gluon_a.dot(pb)

A_s = (term_1 + 2.0 * term_2 - 2.0 * term_3) / pys.s(pa, pb)

print "   -> A_s = ", A_s[pol]

#
# t-channel
#

term_4  = pys.Current(u_2.bar(), 'mu', u_a('+')).dot(Gluon_b) * pys.Current(u_a.bar('+'), 'mu', v_1).dot(Gluon_a)
term_4 += pys.Current(u_2.bar(), 'mu', u_a('-')).dot(Gluon_b) * pys.Current(u_a.bar('-'), 'mu', v_1).dot(Gluon_a)

A_t = - pys.i * term_4 / pys.s(pa, p1)

print "   -> A_t = ", A_t[pol]

#
# u-channel
#

term_5  = pys.Current(u_2.bar(), 'mu', u_a('+')).dot(Gluon_a) * pys.Current(u_a.bar('+'), 'mu', v_1).dot(Gluon_b)
term_5 += pys.Current(u_2.bar(), 'mu', u_a('-')).dot(Gluon_a) * pys.Current(u_a.bar('-'), 'mu', v_1).dot(Gluon_b)
term_6  = pys.Current(u_2.bar(), 'mu', u_2('+')).dot(Gluon_a) * pys.Current(u_2.bar('+'), 'mu', v_1).dot(Gluon_b)
term_6 += pys.Current(u_2.bar(), 'mu', u_2('-')).dot(Gluon_a) * pys.Current(u_2.bar('-'), 'mu', v_1).dot(Gluon_b)

A_u = pys.i * (term_5 - term_6) / pys.s(pa, p2)

print "   -> A_u = ", A_u[pol]

#
# Form colour flow amplitudes for each polarisation
#

j0 = []
j1 = []

for pol in pols:

	j0.append((- pys.i * A_s + A_u)[pol])
	j1.append((  pys.i * A_s + A_t)[pol])

#
# Sum over the squared amplitudes with the colour factors and couplings
#

Sum = 0.0
for i, pol in enumerate(pols):

	print "\n   -> Matrix Element squared for polarisation", pol, " = ",

	temp = pys.g_s ** 4 * (16.0 * pys.abs2(j0[i]) + 16.0 * pys.abs2(j1[i]) - 4.0 * (j0[i] * j1[i].conjugate()).real) / 3.0

	# Priunt indiviual amplitudes
	print temp

	# Add the contribution to the total sum
	Sum += temp

print "\n   -> Sum over all polarisation(s) = ",  Sum, "\n"


