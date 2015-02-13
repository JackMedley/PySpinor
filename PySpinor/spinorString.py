

class spinorString(object):

	def __init__(self, *args):

		# Save the arguments
		self.SString = args

	def __str__(self, latex=False):

		retVal = ''

		for arg in self.SString:

			if isinstance(arg, Spinor):

				if arg.barred:
					retVal += '< i |'
				else:
					retVal += '| i >'

			elif isinstance(arg, Spinor.spinor):

				if arg.barred:
					if arg.pol == '+':
						retVal += '< +i |'
					else:
						retVal += '< -i |'
				else:
					if arg.pol == '+':
						retVal += '| i+ >'
					else:
						retVal += '| i- >'

			elif isinstance(arg, str):

				retVal += 'g[' + arg + ']'

			else:
				print '?!'

		return retVal