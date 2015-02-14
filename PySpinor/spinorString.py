from collections import Counter
from Spinor import Spinor
import os.path
import subprocess

class spinorString(object):

	# Pass with dictionary of terms *in order* and a numerical prefactor
	def __init__(self, *args, **kwargs):

		# Empty string where the latex output will be stored
		self.output = ''

		if 'colour' not in kwargs:
			self.SU3 = False
		else:
			self.SU3 = True

		# Save the arguments
		self.String    = args
		self.prefactor = kwargs['prefactor']
		if self.SU3:
			self.colour    = kwargs['colour']

	def clearString(self):

		self.output = ''

	def addString(self):

		# If we already have a string to which we are adding
		output = '$' + self.output.translate(None, '$')

		if self.SU3:
			output += self.prefactor + '[' + self.colour + ']\left('
		else:
			output += self.prefactor + '\left('

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
					output += '\gamma^{\\' + arg + '}'
				else:
					output += '\gamma_{\\' + arg + '}'

				indices.append(arg)

			else:
				print '?!'

		output += '\\right)$'

		# Save and return this string
		self.output = output

		return output

	def getLaTeX(self, fileName):

		self.fileName = fileName

		if os.path.isfile(fileName):
			subprocess.call(["rm", os.path.dirname(os.path.realpath(__file__)) + "/" + fileName + ".tex"])

		# Make the output file
		fileO = open(os.path.dirname(os.path.realpath(__file__)) + "/" + fileName + ".tex", "w")

		fileO.write("""
		\documentclass{article}
		\usepackage{graphicx}
		\\begin{document}\n\n""" + self.output + """\n
		\end{document}""")

	def compileLaTeX(self):

		subprocess.call(["pdflatex", os.path.dirname(os.path.realpath(__file__)) + "/" + self.fileName  + ".tex",  "> /dev/null"])

		subprocess.call(["rm", os.path.dirname(os.path.realpath(__file__)) + "/" + self.fileName + ".aux ",
			               os.path.dirname(os.path.realpath(__file__)) + "/" + self.fileName + ".log"])


