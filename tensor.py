#! /usr/bin/env python

# Import spinor library
from PySpinor2 import *

#
# Define momenta in the problem
#

pa = Momenta(7.500000e+02,  0.000000e+00,  0.000000e+00,  7.500000e+02, incoming=True)  # Incoming gluon
pb = Momenta(7.500000e+02,  0.000000e+00,  0.000000e+00, -7.500000e+02, incoming=True)  # Incoming gluon
p1 = Momenta(7.500000e+02,  1.618995e+02,  6.492524e+02, -2.912573e+02)		        # Outgoing tx
p2 = Momenta(7.500000e+02, -1.618995e+02, -6.492524e+02,  2.912573e+02)		        # Outgoing t

# Check the phase space point conserves momentum
assert Momenta.Conserved() is True, "ERROR! Momentum not conserved..."

# Now things are more complicated since the tops have a non-negligible mass
# Start by decomposing the top momenta into two basis momenta (each)
k1, k2 = p1.decompose()
k3, k4 = p2.decompose()

#
# Define physical spinors
#

# Incoming spinors
u_a = Spinor(pa)
u_b = Spinor(pb)

# Outgoing spinors - there are massive spinors which are constructed by decomposing the
# massive momentum, using the massless spinor conventions and then combined.
u_1 = Spinor(k1)
u_2 = Spinor(k2)
u_3 = Spinor(k3)
u_4 = Spinor(k4)

# Define the gluon currents using these reference momenta
Gluon_a = Gluon(pa, k2, 'nu')['+']
Gluon_b = Gluon(pb, k3, 'mu')['-']

temp = Current(u_a.bar(), 'mu', u_b)

# Test out tensor stuff
print "1st Term\n", Tensor(u_1.bar('+'), 'mu', 'nu', u_2('-')).contract(p1('mu'), p2('nu'))
print "2nd Term\n", Tensor(u_1.bar('+'), 'mu', 'nu', u_2('-')).contract(p2('nu'), p1('mu'))
print "3nd Term\n", Tensor(u_1.bar('+'), 'mu', 'nu', u_2('-')).contract(p1('nu'), p1('mu')) / ((u_1 // u_2) * p1.mass2())

print "4th Term\n", Tensor(u_1.bar('+'), 'mu', 'nu', u_2('-')).contract(p1('mu'), p2('nu')) + \
                    Tensor(u_1.bar('+'), 'mu', 'nu', u_2('-')).contract(p1('nu'), p2('mu'))
print "4th Term\n", Tensor(u_1.bar('+'), 'mu', 'nu', u_2('-')).contract(p1('mu'), p2('nu')) + \
                    Tensor(u_1.bar('+'), 'mu', 'nu', u_2('-')).contract(p2('mu'), p1('nu'))
print "should be...\n", s(p1, p2) * (u_1 // u_2)

term_5 = Tensor(u_4.bar('-'), 'mu', 'nu', u_2('+')).contract(Gluon_b('mu'), Gluon_a('nu'))
# print "term_5", term_5
# term_5 = Tensor(u_4.bar('-'), 'mu', 'nu', u_2('+')).contract(Gluon_a('nu'), Gluon_b('mu'))
# print "term_5", term_5
# term_6 = Tensor(u_3.bar('+'), 'mu', 'nu', u_1('-')).contract(Gluon_b('mu'), Gluon_a('nu'))
# print "term_6", term_6
# term_6 = Tensor(u_3.bar('+'), 'mu', 'nu', u_1('-')).contract(Gluon_a('nu'), Gluon_b('mu'))
# print "term_6", term_6

print "here"
term_6 = Tensor(u_4.bar('-'), 'mu', 'nu', u_2('+')).contract(temp('mu'), Gluon_a('nu'))
