#!/bin/python

# Import spinor library
from PySpinor import *

#
# First MadGraph Test!
# Comparing my calculation of the polarisation summed/averaged ME for e+ e- scattering
# to that obtained from MadGraphs standalone (c++) result.  The exact result is
#
# Amp = e^4 / q^4 * < 4 | mu | 3 >< 1 | mu | 2 >
#
# and for the phase space point below MadGraph returns:
# 'Matrix element = 1.042635e-02 GeV^0'
#
# JM
#

# Define momenta in the problem
p1 = Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00, 7.500000e+02, incoming=True)  # Incoming positron
p2 = Momenta(7.500000e+02, 0.000000e+00, 0.000000e+00,-7.500000e+02, incoming=True)  # Incoming electron
p3 = Momenta(7.500000e+02, 1.663864e+02, 6.672462e+02,-2.993294e+02)		     # Outgoing positron
p4 = Momenta(7.500000e+02,-1.663864e+02,-6.672462e+02, 2.993294e+02)		     # Outgoing electron

# Define spinors
v_1 = Spinor(p1, fermion=False)
u_2 = Spinor(p2)
v_3 = Spinor(p3, fermion=False)
u_4 = Spinor(p4)

# Momentum transfer vector
q = p1 + p2

# Value of the prefactors
prefactor = e ** 4 / q.dot(q) ** 2

# Calculate all combinations of polarisations
cur1 = Current(u_4.bar(), "mu", v_3).dot(Current(v_1.bar(), "mu", u_2))

# Perform all of the contractions
contraction  = abs2(cur1)

# Combine prefactors and the sum of helicity combinations and then average over spins
Amplitude = (prefactor * contraction.sum()) / 4.0

# Return the result
print Amplitude, "GeV"
