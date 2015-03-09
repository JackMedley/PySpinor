from Common import *
from Momenta import Momenta
from Spinor  import Spinor
from Current import Current

## Gluon class
class Gluon(object):

	def __init__(self, p, r, index, upper=True):


		assert isinstance(p, Momenta) is True and \
		       isinstance(r, Momenta) is True and \
		       isinstance(index, str) is True, "ERROR! Wrong Gluon init. arguments..."

		assert p.isMassless() and r.isMassless(), "ERROR! Cant make gluons with massive reference momenta..."

		assert p.dot(r) != 0.0, "ERROR! Physical and reference momenta are collinear..."

		self.p     = p
		self.r     = r
		self.index = index
		self.upper = upper

		self.vector = {}

		u_p = Spinor(p)
		u_r = Spinor(r)

		self.vector['+'] = - Current.current(u_p.bar('-'), index, u_r('-')) / (sqrt(2) * (u_r // u_p))
		self.vector['-'] =   Current.current(u_r.bar('-'), index, u_p('-')) / (sqrt(2) * (u_r ** u_p))

	def __str__(self):

		return '+:\n' + self.vector['+'].__str__() + '\n-:\n' + self.vector['-'].__str__()

	def __getitem__(self, pol):

		assert pol == '+' or pol == '-', "ERROR! Wrong gluon polarisation..."

		if pol == '+':
			return self.vector['+']
		elif pol == '-':
			return self.vector['-']