#! /usr/bin/env python

from PySpinor import *

# Momenta
pa = Momenta(7.500000e+02,  0.000000e+00,  0.000000e+00,  7.500000e+02, incoming=True)
pb = Momenta(7.500000e+02,  0.000000e+00,  0.000000e+00, -7.500000e+02, incoming=True)
p1 = Momenta(7.500000e+02,  1.618995e+02,  6.492524e+02, -2.912573e+02)
p2 = Momenta(7.500000e+02, -1.618995e+02, -6.492524e+02,  2.912573e+02)

# Spinors
u_a = Spinor(pa)
u_b = Spinor(pb)
u_1 = Spinor(p1)
u_2 = Spinor(p2)

# Trial spinor string
x = spinorString(u_a.bar('+'), 'mu', u_b)
print x