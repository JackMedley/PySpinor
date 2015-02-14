from collections import Counter
from Spinor import Spinor

class spinorString(object):

	# Pass with dictionary of terms *in order* and a numerical prefactor
	def __init__(self, *args, **kwargs):

		if 'colour' not in kwargs:
			self.SU3 = False
		else:
			self.SU3 = True

		# Save the arguments
		self.String    = args
		self.prefactor = kwargs['prefactor']
		if self.SU3:
			self.colour    = kwargs['colour']

	def getString(self):

		if self.SU3:
			output = '$' + self.prefactor + '[' + self.colour + ']\left('
		else:
			output = '$' + self.prefactor + '\left('


		indices = []

		for arg in self.String:

			if isinstance(arg, Spinor):

				if arg.fermion:
					if arg.barred:
						output += '\overline{u}(p_{' + str(arg.momentum.label) + '})'
					else:
						output += 'u(p_{' + str(arg.momentum.label) + '})'
				else:
					if arg.barred:
							output += '\overline{v}(p_{' + str(arg.momentum.label) + '})'
					else:
						output += 'v(p_{' + str(arg.momentum.label) + '})'

			elif isinstance(arg, Spinor.spinor):

				if arg.fermion:
					if arg.barred:
						if arg.pol == '+':
							output += '\overline{u}^+(p_{' + str(arg.momentum.label) + '})'
						else:
							output += '\overline{u}^-(p_{' + str(arg.momentum.label) + '})'
					else:
						if arg.pol == '+':
							output += 'u^+(p_{' + str(arg.momentum.label) + '})'
						else:
							output += 'u^-(p_{' + str(arg.momentum.label) + '})'
				else:
					if arg.barred:
						if arg.pol == '+':
							output += '\overline{v}^+(p_{' + str(arg.momentum.label) + '})'
						else:
							output += '\overline{v}^-(p_{' + str(arg.momentum.label) + '})'
					else:
						if arg.pol == '+':
							output += 'v^+(p_{' + str(arg.momentum.label) + '})'
						else:
							output += 'v^-(p_{' + str(arg.momentum.label) + '})'

			elif isinstance(arg, str):

				assert Counter(indices)[arg] < 2, "ERROR! Index " + arg + " appears more than once in spinor string"

				if arg not in indices:
					output += '\gamma^{' + arg + '}'
				else:
					output += '\gamma_{' + arg + '}'

				indices.append(arg)

			else:
				print '?!'

		output += '$'

		return output