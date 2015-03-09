from Common import *
from Momenta import Momenta
from LorentzVector import LorentzVector

## A functor class whick holds an item (a current or a LorentzVector, ...) for each polarisation in the problem
class polDictionary(object):

	# Flag to avoid multiple calls to 'getPolarisations'
	Init     = False
	baseDict = None

	def __init__(self, Dict):

		# Save the dictionary
		self.polarisations = Dict

		# Run the zero checker on it
		self.zeroChecker()

	def real(self):

		copy = deepcopy(self)

		for key in copy.polarisations:

			copy.polarisations[key] = copy.polarisations[key].real

		return copy

	def imag(self):

		copy = deepcopy(self)

		for key in copy.polarisations:

			copy.polarisations[key] = copy.polarisations[key].imag

		return copy

	## Zero(ish) checker to fix screwed up numerics where things are basically zero but not quite
	def zeroChecker(self):

		if isinstance(self.polarisations.itervalues().next(), float) or \
		   isinstance(self.polarisations.itervalues().next(), int)   or \
		   isinstance(self.polarisations.itervalues().next(), complex):

		   	# Loop over dictionary
			for key in self.polarisations:

				if abs(self.polarisations[key]) < TOLERANCE:

					self.polarisations[key] = zero

		if isinstance(self.polarisations.itervalues().next(), LorentzVector):

		   	# Loop over dictionary
			for key in self.polarisations:

				if abs(self.polarisations[key].T()) < TOLERANCE:

					self.polarisations[key].T(zero)

				if abs(self.polarisations[key].X()) < TOLERANCE:

					self.polarisations[key].X(zero)

				if abs(self.polarisations[key].Y()) < TOLERANCE:

					self.polarisations[key].T(zero)

				if abs(self.polarisations[key].Z()) < TOLERANCE:

					self.polarisations[key].Z(zero)

	## Define functor behaviour
	def __call__(self, index):

		copy = deepcopy(self)

		for key in copy.polarisations:

			if isinstance(copy.polarisations[key], LorentzVector) == LorentzVector:

				copy.polarisations[key](index)

		return copy

	def __getitem__(self, i):

		return self.polarisations[i]

	def __neg__ (self):

		copy = deepcopy(self)

		for key in copy.polarisations:

			copy.polarisations[key] *= -1.0

		return copy

	def __pos__ (self):

		return self

	def __iadd__(self, other):

		assert isinstance(other, polDictionary), "ERROR!  Cant add non-polDictionary to polDictionary..."

		temp = {}

		for key in self.polarisations:

			temp[key] = self.polarisations[key] + other.polarisations[key]

		return polDictionary(temp)

	def __str__(self):

		temp = ''

		for key in self.polarisations:
			temp += str(key) + ' ' + self.polarisations[key].__str__() + '\n'

		return temp

	def __pow__(self, other):

		temp = {}

		for key in self.polarisations:

			temp[key] = self.polarisations[key] ** other

		return polDictionary(temp)

	def __eq__(self, other):

		Pass = True

		for key in self.polarisations:

			if isinstance(other, float) or isinstance(other, int) or isinstance(other, complex):

				if abs(self.polarisations[key] - other) > TOLERANCE:

					Pass = False

			elif isinstance(other, polDictionary):

				if abs(self.polarisations[key] - other.polarisations[key]) > TOLERANCE:

					Pass = false

		return Pass

	def __add__(self, other):

		temp = {}

		for key in self.polarisations:
			temp[key] = self.polarisations[key] + other.polarisations[key]

		return polDictionary(temp)

	def __sub__(self, other):

		temp = {}

		for key in self.polarisations:
			temp[key] = self.polarisations[key] - other.polarisations[key]

		return polDictionary(temp)

	def __mul__(self, other):

		temp = {}

		for key in self.polarisations:

			if isinstance(other, polDictionary):
				temp[key] = self.polarisations[key] * other.polarisations[key]
			elif isinstance(other, float) or isinstance(other, int) or isinstance(other, complex):
				temp[key] = self.polarisations[key] * other

		return polDictionary(temp)

	def __rmul__(self, other):

		temp = {}

		# Construct the multiplied dictionary
		for key in self.polarisations:

			if isinstance(other, polDictionary):
				temp[key] = other.polarisations[key] * self.polarisations[key]
			if isinstance(other, float) or isinstance(other, int) or isinstance(other, complex):
				temp[key] = self.polarisations[key] * other

		return polDictionary(temp)

	def __div__(self, other):

		assert isinstance(other, polDictionary) is True or \
		       isinstance(other, float)         is True or \
		       isinstance(other, complex)       is True or \
		       isinstance(other, int)           is True, "ERROR!"

		temp = {}

		for key in self.polarisations:

			if isinstance(other, polDictionary):
					temp[key] = self.polarisations[key] / other.polarisations[key]
			else:
				for key in self.polarisations:
					temp[key] = self.polarisations[key] / other

		return polDictionary(temp)

	def dot(self, other):

		temp = {}

		assert isinstance(other, LorentzVector) or \
		       isinstance(other, polDictionary), "ERROR! Wrong args. in dot..." + str(type((other)))

		# For each polarisation
		for key in self.polarisations:


			if isinstance(other, LorentzVector):

				temp[key] = self.polarisations[key].dot(other)

			else:
				temp[key] = self.polarisations[key].dot(other.polarisations[key])

		return polDictionary(temp)

	def conjugate(self):

		temp = {}

		for key in self.polarisations:
			temp[key] = self.polarisations[key].conjugate()
		return polDictionary(temp)

	@classmethod
	def setPolarisations(self, List):

		return

	@classmethod
	def getPolarisations(self, List=False, ForceNew=False, Fixed=None):

		# If we've already generated the dictionary then return a clean copy rather than recalculating
		if polDictionary.Init == True and polDictionary.baseDict != None and ForceNew == False:

			if List == False:
				return copy.deepcopy(polDictionary.baseDict)

			else:
				temp = []

				for key in polDictionary.baseDict:

					temp.append(key)

				return temp

		# Make an empty dictionary
		temp = {}

		# Generate all the polarisations
		keys = list(product(('-', '+'), repeat=Momenta.pCount))

		# Loop over the list and add each one as a key
		for key in keys:
			temp[key] = None

		# Set the flag to say we've run this once
		polDictionary.Init = True

		# Save the dictionary
		polDictionary.baseDict = temp

		if List:

			# Return a list of polarisations
			temp = []

			for key in polDictionary.baseDict:

				temp.append(key)

			return temp
		else:
			# Return the dictionary
			return temp