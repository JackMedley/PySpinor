from Common import *
from polDictionary import polDictionary
from LorentzVector import LorentzVector


## Class storing spinor current information
class Current(polDictionary):

	@classmethod
	def initFromDict(cls, instance, Dict):

		temp = deepcopy(instance)
		temp.polarisations = Dict
		return temp

	## Current constructor
	def __init__(self, spin1, index, spin2, upper=True, Dict=None):
		from Spinor import Spinor

		assert isinstance(spin1, Spinor)        is True or  \
		       isinstance(spin1, Spinor.spinor) is True and \
		       isinstance(spin2, Spinor)        is True or  \
		       isinstance(spin2, Spinor.spinor) is True, "ERROR! Wrong argument for Current..."

		assert isinstance(index, str) is True, "ERROR! Wrong index argument for Current..." + str(type(index))

		self.index = index
		self.upper = upper
		self.spin1 = spin1
		self.spin2 = spin2

		# Get a dictionary of all possible polarisation combinations
		if Dict != None:
			self.polarisations = Dict
			polDictionary.__init__(self, self.polarisations)

		# Unless we are init-ing from a dictionary
		else:
			self.polarisations = polDictionary.getPolarisations()

			# Fill in all the possible polarisations depending on the arguments
			if isinstance(spin1, Spinor) and isinstance(spin2, Spinor):

				for key in self.polarisations:

					if   key[spin1.ID] == '+' and key[spin2.ID] == '+':
						self.polarisations[key] = Current.current(spin1('+'), index, spin2('+'), upper=upper)
					elif key[spin1.ID] == '-' and key[spin2.ID] == '+':
						self.polarisations[key] = Current.current(spin1('-'), index, spin2('+'), upper=upper)
					elif key[spin1.ID] == '+' and key[spin2.ID] == '-':
						self.polarisations[key] = Current.current(spin1('+'), index, spin2('-'), upper=upper)
					elif key[spin1.ID] == '-' and key[spin2.ID] == '-':
						self.polarisations[key] = Current.current(spin1('-'), index, spin2('-'), upper=upper)
					else:
						print "CHECK THIS"

			# If spin1 is a reference spinor
			elif isinstance(spin1, Spinor.spinor) and isinstance(spin2, Spinor):

				# Loop over all polarisations
				for key in self.polarisations:

					# Calculate for all the possible spin2 polarisations while being
					# careful NOT to enforce spin1's physical polarisation since it
					# could be that we are summing over an unphysical spinor
					if key[spin2.ID] == '+':
						self.polarisations[key] = Current.current(spin1, index, spin2('+'), upper=upper)
					elif key[spin2.ID] == '-':
						self.polarisations[key] = Current.current(spin1, index, spin2('-'), upper=upper)
					else:
						print "CHECK THIS"

			# If spin1 is a reference spinor
			elif isinstance(spin1, Spinor) and isinstance(spin2, Spinor.spinor):

				# Loop over all polarisations
				for key in self.polarisations:

					# Calculate for all the possible spin1 polarisations while being
					# careful NOT to enforce spin2's physical polarisation since it
					# could be that we are summing over an unphysical spinor
					if key[spin1.ID] == '+':
						self.polarisations[key] = Current.current(spin1('+'), index, spin2, upper=upper)
					elif key[spin1.ID] == '-':
						self.polarisations[key] = Current.current(spin1('-'), index, spin2, upper=upper)
					else:
						print "CHECK THIS"

			# Else they are both reference spinors
			else:
				for key in self.polarisations:
					self.polarisations[key] = Current.current(spin1, index, spin2, upper=upper)

		# pass the object to the parent constructor
		polDictionary.__init__(self, self.polarisations)

	def __rmul__(self, other):

		temp = {}

		# Construct the multiplied dictionary
		for key in self.polarisations:

			if isinstance(other, polDictionary):
				temp[key] = other.polarisations[key] * self.polarisations[key]
			if isinstance(other, float) or isinstance(other, int) or isinstance(other, complex):
				temp[key] = self.polarisations[key] * other

		return Current(self.spin1, self.index, self.spin2, upper=self.upper, Dict=temp)

	def __div__(self, other):

		temp = {}

		# Construct the divided dictionary
		for key in self.polarisations:

			if isinstance(other, polDictionary):
				temp[key] = other.polarisations[key] / self.polarisations[key]
			if isinstance(other, float) or isinstance(other, int) or isinstance(other, complex):
				temp[key] = self.polarisations[key] / other

		return Current(self.spin1, self.index, self.spin2, upper=self.upper, Dict=temp)

	## Lorentz vector currents
	class current(LorentzVector):

		# Keep track of all the indices in the calculatation
		indices = {}

		def __init__ (self, spin1, index, spin2, upper=True, Vector=None):
			from Spinor import Spinor

			assert isinstance(spin1, Spinor.spinor) is True and \
			       isinstance(spin1, Spinor.spinor) is True, "ERROR! Wrong argument for current..."

			self.spin1 = spin1
			self.spin2 = spin2
			self.upper = upper
			self.index = index

			if Vector != None:
				vector = Vector
			else:
				vector = np.matrix([(spin1 * gamma[0] * spin2).item(0),
				                    (spin1 * gamma[1] * spin2).item(0),
				                    (spin1 * gamma[2] * spin2).item(0),
				                    (spin1 * gamma[3] * spin2).item(0)])

			LorentzVector.__init__(self, vector.item(0), vector.item(1), vector.item(2), vector.item(3), upper=True, index=self.index)