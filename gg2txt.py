#! /usr/bin/env python

# Import spinor library
from PySpinor import *

#
# Define momenta in the problem
#

pa = Momenta(7.500000e+02,  0.000000e+00,  0.000000e+00,  7.500000e+02, incoming=True)  # Incoming gluon
pb = Momenta(7.500000e+02,  0.000000e+00,  0.000000e+00, -7.500000e+02, incoming=True)  # Incoming gluon
p1 = Momenta(7.500000e+02,  1.618995e+02,  6.492524e+02, -2.912573e+02)		              # Outgoing tx
p2 = Momenta(7.500000e+02, -1.618995e+02, -6.492524e+02,  2.912573e+02)		              # Outgoing t

# Check the phase space point conserves momentum
assert Momenta.Conserved() is True, "ERROR! Momentum not conserved..."

# Now things are more complicated since the tops have a non-negligible mass
# Start by decomposing the top momenta into two basis momenta (each)
k1, k2 = p1.decompose()
k3, k4 = p2.decompose()

# Spinors
u_a = Spinor(pa)
u_b = Spinor(pb)
u_1 = Spinor(k1)
u_2 = Spinor(k2)
u_3 = Spinor(k3)
u_4 = Spinor(k4)

polarisations = [('+', '-', '+', '-')]
pol = polarisations[0]

# Choose two massless reference momenta - changes in this should be equivalent to a gauge transformation
ra = k2
rb = k3

# Define the gluon currents using these reference momenta - functor behaviour defines polarisation
Gluon_a = Gluon(pa, ra, 'nu')['+']
Gluon_b = Gluon(pb, rb, 'mu')['-']

# Useful constants
f1 = (u_3 // u_4) / p1.mass()
f2 = (u_2 // u_1) / p1.mass()
m  = p1.mass()

#
# s-channel
#

# Calculation using massless spinors
term_1 =       f1 * Current(u_4.bar('-'), 'mu', u_1('-')).dot(pb - pa) * Gluon_a.dot(Gluon_b)
term_2 =       f2 * Current(u_2.bar('-'), 'mu', u_3('-')).dot(pb - pa) * Gluon_a.dot(Gluon_b)
term_3 = 2.0 * f1 * Current(u_4.bar('-'), 'mu', u_1('-')).dot(Gluon_a) * Gluon_b.dot(pa)
term_4 = 2.0 * f2 * Current(u_2.bar('-'), 'mu', u_3('-')).dot(Gluon_a) * Gluon_b.dot(pa)
term_5 = 2.0 * f1 * Current(u_4.bar('-'), 'mu', u_1('-')).dot(Gluon_b) * Gluon_a.dot(pb)
term_6 = 2.0 * f2 * Current(u_2.bar('-'), 'mu', u_3('-')).dot(Gluon_b) * Gluon_a.dot(pb)

# Combine terms
A_s =  (- term_1 + term_2 - term_3 + term_4 + term_5 - term_6) / s(pa, pb)

print "\n   -> A_s   = ", A_s[pol]

#
# t-channel
#

term_1 = f1 * Current(u_4.bar('-'), 'mu', u_a('-')).dot(Gluon_b) * \
              Current(u_a.bar('-'), 'mu', u_1('-')).dot(Gluon_a)
term_2 = f1 * Current(u_4.bar('-'), 'mu', u_1('-')).dot(Gluon_b) * \
              Current(u_1.bar('-'), 'mu', u_1('-')).dot(Gluon_a)
term_3 = f1 * Current(u_4.bar('-'), 'mu', u_2('-')).dot(Gluon_b) * \
              Current(u_2.bar('-'), 'mu', u_1('-')).dot(Gluon_a)
term_4 = f2 * Current(u_a.bar('-'), 'mu', u_3('-')).dot(Gluon_b) * \
              Current(u_2.bar('-'), 'mu', u_a('-')).dot(Gluon_a)
term_5 = f2 * Current(u_1.bar('-'), 'mu', u_3('-')).dot(Gluon_b) * \
              Current(u_2.bar('-'), 'mu', u_1('-')).dot(Gluon_a)
term_6 = f2 * Current(u_2.bar('-'), 'mu', u_3('-')).dot(Gluon_b) * \
              Current(u_2.bar('-'), 'mu', u_2('-')).dot(Gluon_a)

term_7 = f1 * f2 * m * Tensor(u_4.bar('-'), 'mu', 'nu', u_2('+')).contract(Gluon_b('mu'), Gluon_a('nu'))
term_8 =           m * Tensor(u_3.bar('+'), 'mu', 'nu', u_1('-')).contract(Gluon_b('mu'), Gluon_a('nu'))

A_t = - i * (-term_1 + term_2 + term_3 + term_4 - term_5 - term_6 + term_7 - term_8) / s(pa, p1)

print "   -> A_t   = ", A_t[pol]

#
# u-channel
#

term_1 = f1 * Current(u_4.bar('-'), 'mu', u_a('-')).dot(Gluon_a) * \
              Current(u_a.bar('-'), 'mu', u_1('-')).dot(Gluon_b)
term_2 = f1 * Current(u_4.bar('-'), 'mu', u_3('-')).dot(Gluon_a) * \
              Current(u_3.bar('-'), 'mu', u_1('-')).dot(Gluon_b)
term_3 = f1 * Current(u_4.bar('-'), 'mu', u_4('-')).dot(Gluon_a) * \
              Current(u_4.bar('-'), 'mu', u_1('-')).dot(Gluon_b)
term_4 = f2 * Current(u_a.bar('-'), 'mu', u_3('-')).dot(Gluon_a) * \
              Current(u_2.bar('-'), 'mu', u_a('-')).dot(Gluon_b)
term_5 = f2 * Current(u_3.bar('-'), 'mu', u_3('-')).dot(Gluon_a) * \
              Current(u_2.bar('-'), 'mu', u_3('-')).dot(Gluon_b)
term_6 = f2 * Current(u_4.bar('-'), 'mu', u_3('-')).dot(Gluon_a) * \
              Current(u_2.bar('-'), 'mu', u_4('-')).dot(Gluon_b)
term_7 = m * f1 * f2 * Tensor(u_4.bar('-'), 'nu', 'mu', u_2('+')).contract(Gluon_b('mu'), Gluon_a('nu'))
term_8 = m           * Tensor(u_3.bar('+'), 'nu', 'mu', u_1('-')).contract(Gluon_b('mu'), Gluon_a('nu'))

A_u =  i * (-term_1 + term_2 + term_3 + term_4 - term_5 - term_6 + term_7 - term_8) / s(pa, p2)

print "   -> A_u   = ", A_u[pol]

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

	# Priunt individual amplitudes
	print temp

	# Add the contribution to the total sum
	itr += 1
	Sum += temp

print "\n"
if itr > 1:
	print "   -> Sum over all polarisation(s) = ",  Sum, "\n"
