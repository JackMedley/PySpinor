from PySpinor import *

elements = np.array([[1.0,  2.0,  3.0,  4.0],  \
	                 [5.0,  6.0,  7.0,  8.0],  \
	                 [9.0,  10.0, 11.0, 12.0], \
	                 [13.0, 14.0, 15.0, 16.0] ])
indices = ['a', 'b']
ulList = ['upper', 'upper']

T = Tensor(elements, indices, ulList)

# print -T+T
# print T.GetIndices()

p = Tensor(np.array([1.0, 2.0, 3.0, 4.0]), ['b'], ['upper'])

pa = Tensor(np.array([1.0, 2.0, 3.0, 4.0]), ['a'], ['upper'])

p.LowerIndex('b')

print p

print T*p
