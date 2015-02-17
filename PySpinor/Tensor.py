from Common import *
from polDictionary import polDictionary
from LorentzVector import LorentzVector
from Spinor import Spinor

class Tensor(polDictionary):

	def __init__(self, spin1, index1, index2, spin2, upper1=True, upper2=True):

		assert isinstance(spin1, Spinor)        is True or  \
		       isinstance(spin1, Spinor.spinor) is True and \
		       isinstance(spin2, Spinor)        is True or  \
		       isinstance(spin2, Spinor.spinor) is True, "ERROR! Wrong argument for Current..."

		self.spin1  = spin1
		self.spin2  = spin2
		self.index1 = index1
		self.index2 = index2
		self.upper1 = upper1
		self.upper2 = upper2

		# Generate a clean dictionary of all possible polarisation combinations
		self.polarisations = polDictionary.getPolarisations()

		# Fill in all the possible polarisations depending on the arguments
		if isinstance(spin1, Spinor) and isinstance(spin2, Spinor):

			for key in self.polarisations:

				if key[spin1.ID] == '+' and key[spin2.ID] == '+':
					self.polarisations[key] = Tensor.tensor(spin1('+'), index1, index2, spin2('+'))
				elif key[spin1.ID] == '+' and key[spin2.ID] == '-':
					self.polarisations[key] = Tensor.tensor(spin1('+'), index1, index2, spin2('-'))
				elif key[spin1.ID] == '-' and key[spin2.ID] == '+':
					self.polarisations[key] = Tensor.tensor(spin1('-'), index1, index2, spin2('+'))
				elif key[spin1.ID] == '-' and key[spin2.ID] == '-':
					self.polarisations[key] = Tensor.tensor(spin1('-'), index1, index2, spin2('-'))
				else:
					print "CHECK THIS"

		# If spin1 is a reference spinor
		elif isinstance(spin1, Spinor.spinor) and isinstance(spin2, Spinor):

			# Loop over all polarisations
			for key in self.polarisations:

				if key[spin2.ID] == '+':
					self.polarisations[key] = Tensor.tensor(spin1, index1, index2, spin2('+'))
				elif key[spin2.ID] == '-':
					self.polarisations[key] = Tensor.tensor(spin1, index1, index2, spin2('-'))
				else:
					print "CHECK THIS"

		# If spin1 is a reference spinor
		elif isinstance(spin1, Spinor) and isinstance(spin2, Spinor.spinor):

			# Loop over all polarisations
			for key in self.polarisations:

				if key[spin1.ID] == '+':
					self.polarisations[key] = Tensor.tensor(spin1('+'), index1, index2, spin2)

				elif key[spin1.ID] == '-':
					self.polarisations[key] = Tensor.tensor(spin1('-'), index1, index2, spin2)
				else:
					print "CHECK THIS"

		# Else they are both reference spinors
		else:
			for key in self.polarisations:
				self.polarisations[key] = Tensor.tensor(spin1, index1, index2, spin2)

	# Contraction for the full polDictionary Tensor object
	def contract(self, other1, other2):

		assert isinstance(other1, polDictionary) or isinstance(other1, LorentzVector) and \
		       isinstance(other2, polDictionary) or isinstance(other2, LorentzVector), "Cant contract Tensor with non-polDictionary type..."

		# Index checks
		assert other1.index != other2.index, "Wrong indices for tensory contraction..."
		assert other1.index == self.index1 or \
		       other1.index == self.index2,  "Wrong indices for tensory contraction..."
		assert other2.index == self.index1 or \
		       other2.index == self.index2,  "Wrong indices for tensory contraction..."

		temp = {}

		for key in self.polarisations:

			tempVar = 0.0
			sign    = +1

			# Loop over tensor
			for i in range(0, 4):

				for j in range(0, 4):

					if i == j == 0:
						sign = +1.0
					elif i == 0 or j == 0:
						sign = -1.0
					else:
						sign = +1.0

					if other1.index == self.index1:
						if isinstance(other1, polDictionary) and isinstance(other2, polDictionary):
							tempVar += sign * self[key].index1Vec[i][j] * other1[key][j] * other2[key][i]
						elif isinstance(other1, polDictionary):
							tempVar += sign * self[key].index1Vec[i][j] * other1[key][j] * other2[i]
						elif isinstance(other2, polDictionary):
							tempVar += sign * self[key].index1Vec[i][j] * other1[j] * other2[key][i]
						else:
							tempVar += sign * self[key].index1Vec[i][j] * other1[j] * other2[i]

					elif other2.index == self.index1:
						if isinstance(other1, polDictionary) and isinstance(other2, polDictionary):
							tempVar += sign * self[key].index1Vec[i][j] * other1[key][i] * other2[key][j]
						elif isinstance(other1, polDictionary):
							tempVar += sign * self[key].index1Vec[i][j] * other1[key][i] * other2[j]
						elif isinstance(other2, polDictionary):
							tempVar += sign * self[key].index1Vec[i][j] * other1[i] * other2[key][j]
						else:
							tempVar += sign * self[key].index1Vec[i][j] * other1[i] * other2[j]

			temp[key] = tempVar

		return polDictionary(temp)

	class tensor:

		def __init__ (self, spin1, index1, index2, spin2, upper1=True, upper2=True):

			assert isinstance(spin1, Spinor.spinor) is True and isinstance(spin1, Spinor.spinor) is True, "ERROR! Wrong argument for current..."

			# Save the index
			self.index1 = index1
			self.index2 = index2

			self.upper1 = upper1
			self.upper2 = upper2

			self.index1Vec = []
			self.index2Vec = []

			# Calculate 'row' formation of tensor - traversed by first index
			self.index1Vec.append(LorentzVector((spin1 * gamma[0] * gamma[0] * spin2).item(0), \
			                                    (spin1 * gamma[0] * gamma[1] * spin2).item(0), \
			                                    (spin1 * gamma[0] * gamma[2] * spin2).item(0), \
			                                    (spin1 * gamma[0] * gamma[3] * spin2).item(0), True, index=index1))

			self.index1Vec.append(LorentzVector((spin1 * gamma[1] * gamma[0] * spin2).item(0), \
			                                    (spin1 * gamma[1] * gamma[1] * spin2).item(0), \
			                                    (spin1 * gamma[1] * gamma[2] * spin2).item(0), \
			                                    (spin1 * gamma[1] * gamma[3] * spin2).item(0), True, index=index1))

			self.index1Vec.append(LorentzVector((spin1 * gamma[2] * gamma[0] * spin2).item(0), \
			                                    (spin1 * gamma[2] * gamma[1] * spin2).item(0), \
			                                    (spin1 * gamma[2] * gamma[2] * spin2).item(0), \
			                                    (spin1 * gamma[2] * gamma[3] * spin2).item(0), True, index=index1))

			self.index1Vec.append(LorentzVector((spin1 * gamma[3] * gamma[0] * spin2).item(0), \
			                                    (spin1 * gamma[3] * gamma[1] * spin2).item(0), \
			                                    (spin1 * gamma[3] * gamma[2] * spin2).item(0), \
			                                    (spin1 * gamma[3] * gamma[3] * spin2).item(0), True, index=index1))

			# Calculate 'column' formation of tensor - traversed by second index
			self.index2Vec.append(LorentzVector((spin1 * gamma[0] * gamma[0] * spin2).item(0), \
			                                    (spin1 * gamma[1] * gamma[0] * spin2).item(0), \
			                                    (spin1 * gamma[2] * gamma[0] * spin2).item(0), \
			                                    (spin1 * gamma[3] * gamma[0] * spin2).item(0), True, index=index2))

			self.index2Vec.append(LorentzVector((spin1 * gamma[0] * gamma[1] * spin2).item(0), \
			                                    (spin1 * gamma[1] * gamma[1] * spin2).item(0), \
			                                    (spin1 * gamma[2] * gamma[1] * spin2).item(0), \
			                                    (spin1 * gamma[3] * gamma[1] * spin2).item(0), True, index=index2))

			self.index2Vec.append(LorentzVector((spin1 * gamma[0] * gamma[2] * spin2).item(0), \
			                                    (spin1 * gamma[1] * gamma[2] * spin2).item(0), \
			                                    (spin1 * gamma[2] * gamma[2] * spin2).item(0), \
			                                    (spin1 * gamma[3] * gamma[2] * spin2).item(0), True, index=index2))

			self.index2Vec.append(LorentzVector((spin1 * gamma[0] * gamma[3] * spin2).item(0), \
			                                    (spin1 * gamma[1] * gamma[3] * spin2).item(0), \
			                                    (spin1 * gamma[2] * gamma[3] * spin2).item(0), \
			                                    (spin1 * gamma[3] * gamma[3] * spin2).item(0), True, index=index2))

		# def __call__(self, id1=None, id2=None):

		# 	if id1 != None and id2 != None:

		# 		if id2 == 0 or id2 == 'T':
		# 			return self.index1Vec[id1].T()
		# 		if id2 == 1 or id2 == 'X':
		# 			return self.index1Vec[id1].X()
		# 		if id2 == 2 or id2 == 'Y':
		# 			return self.index1Vec[id1].Y()
		# 		if id2 == 3 or id2 == 'Z':
		# 			return self.index1Vec[id1].Z()

		# 	elif id2 == None:

		# 		return self.index1Vec[id1]

		# 	elif id2 == None:

		# 		return self.index2Vec[id2]

		# 	else:
		# 		print "Function behaviour requires (at least) one index..."

		# Flip the ordering of the indices
		def transpose(self):

			temp        = self.index1
			self.index1 = self.index2
			self.index2 = temp

			temp           = self.index1Vec
			self.index1Vec = self.index2Vec
			self.index2Vec = temp

			return self

		def __str__(self):

			return "(" + self.index1 + " = 0) " + str(self.index1Vec[0]) + "\n" \
			       "(" + self.index1 + " = 1) " + str(self.index1Vec[1]) + "\n" \
			       "(" + self.index1 + " = 2) " + str(self.index1Vec[2]) + "\n" \
			       "(" + self.index1 + " = 3) " + str(self.index1Vec[3]) + "\n"

		def __neg__ (self):

			temp = deepcopy(self)

			temp.index1Vec[0] *= -1.0
			temp.index1Vec[1] *= -1.0
			temp.index1Vec[2] *= -1.0
			temp.index1Vec[3] *= -1.0
			temp.index2Vec[0] *= -1.0
			temp.index2Vec[1] *= -1.0
			temp.index2Vec[2] *= -1.0
			temp.index2Vec[3] *= -1.0

			return temp