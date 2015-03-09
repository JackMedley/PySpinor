#! /usr/bin/env python

#
# Common file with useful definitions and imports
#

import numpy as np
import copy
from copy import deepcopy
from cmath import *
from itertools import product

# Make sure numpy prints everything correctly
np.set_printoptions(threshold=np.nan)

# Numerical tolerance for calculation
TOLERANCE = 1E-8

# Useful constants and matrices
zero  = complex(0.0, 0.0)
i     = complex(0.0, 1.0)
I4    = np.identity(4)
zero4 = np.matrix([[zero,zero,zero,zero],
	           [zero,zero,zero,zero],
	           [zero,zero,zero,zero],
	           [zero,zero,zero,zero]])

# The electric and strong charges
e   = 3.079538e-01
a_s = 1.180000e-01
g_s = 1.2177157848
pi  = 3.1415926535

# Useful QCD factors
C_f = 4.0 / 3.0
C_a = 3.0

# Define the metric
metric = np.array([[1.0, 0.0, 0.0, 0.0],
	            [0.0,-1.0, 0.0, 0.0],
                    [0.0, 0.0,-1.0, 0.0],
                    [0.0, 0.0, 0.0,-1.0]])

# Define gamma matrices N.B. these have upper indices
g0    = np.matrix([[0.0, 0.0,  1.0,  0.0], [0.0, 0.0,  0.0, 1.0], [1.0,  0.0,  0.0, 0.0], [0.0,  1.0, 0.0,  0.0]])
g1    = np.matrix([[0.0, 0.0,  0.0, -1.0], [0.0, 0.0, -1.0, 0.0], [0.0,  1.0,  0.0, 0.0], [1.0,  0.0, 0.0,  0.0]])
g2    = np.matrix([[0.0, 0.0,  0.0,    i], [0.0, 0.0, -i,   0.0], [0.0, -i,    0.0, 0.0], [i ,   0.0, 0.0,  0.0]])
g3    = np.matrix([[0.0, 0.0, -1.0,  0.0], [0.0, 0.0,  0.0, 1.0], [1.0,  0.0,  0.0, 0.0], [0.0, -1.0, 0.0,  0.0]])
g5    = i * g0 * g1 * g2 * g3
gamma = [g0, g1, g2, g3]

# Define projection operators
w_p = (I4 + g5) / 2.0
w_m = (I4 - g5) / 2.0