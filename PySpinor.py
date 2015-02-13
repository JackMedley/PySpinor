#! /usr/bin/env python

from Common import *
from copy import deepcopy

print "\n"

__array_priority__ = 100

# Functions to find the kinematic invariant s_pq
def s(p, q):

	assert isinstance(p, LorentzVector) is True and \
	       isinstance(q, LorentzVector) is True, "ERROR! s called with massive (or non-) LorentzVector..."

	return 2.0 * p.dot(q)

def abs2(complexVar):

	if isinstance(complexVar, complex):
		return (complexVar * complexVar.conjugate()).real
	else:
		return (complexVar * complexVar.conjugate()).real()

def avgFactor():

	return  2 ** - (Momenta.pCount - Momenta.iCount)

# Parent class for anything with a Lorentz index
class LorentzVector(object):

	def __init__(self, E, X, Y, Z, upper, index = 'mu'):

		# Save the state of the Lorentz index
		self.upper = upper

		# Save the index
		self.index = index

		# Form vector
		self.vector = np.matrix([E, X, Y, Z])

	# Getter Methods
	def T(self, value=None):

		if value == None:
			return self.vector.item(0)
		else:
			self.vector.itemset(0, value)

	def X(self, value=None):

		if value == None:
			return self.vector.item(1)
		else:
			self.vector.itemset(1, value)

	def Y(self, value=None):

		if value == None:
			return self.vector.item(2)
		else:
			self.vector.itemset(2, value)

	def Z(self, value=None):

		if value == None:
			return self.vector.item(3)
		else:
			self.vector.itemset(3, value)

	def conjugate(self):

		return LorentzVector(self.T().conjugate(), self.X().conjugate(), self.Y().conjugate(), self.Z().conjugate(), upper=self.upper, index=self.index)

	# Overloads
	def __call__(self, index):

		assert isinstance(index, str), "LorentzVector can only be called with index..."

		temp = deepcopy(self)
		temp.index = index

		return temp

	def __str__(self):
		return "(" + str(self.T()) + "; " + str(self.X()) + ", " + str(self.Y()) + ", " + str(self.Z()) + ")"

	def __add__(self, other):

		if self.upper != other.upper:
			print "ERROR!  Cant add LorentzVectors with co and contra indices!"

		if isinstance(self, Momenta):
			return Momenta      (self.T() + other.T(), self.X() + other.X(), self.Y() + other.Y(), self.Z() + other.Z(), upper=self.upper, incoming=self.incoming, physical=False)
		else:
			return LorentzVector(self.T() + other.T(), self.X() + other.X(), self.Y() + other.Y(), self.Z() + other.Z(), self.upper)

	def __sub__(self, other):

		if self.upper != other.upper:
			print "ERROR!  Cant subtract LorentzVectors with co and contra indices!"

		if isinstance(self, Momenta):
			return Momenta      (self.T() - other.T(), self.X() - other.X(), self.Y() - other.Y(), self.Z() - other.Z(), upper=self.upper, incoming=self.incoming, physical=False)
		else:
			return LorentzVector(self.T() - other.T(), self.X() - other.X(), self.Y() - other.Y(), self.Z() - other.Z(), self.upper)

	def __neg__ (self):

		temp = deepcopy(self)
		temp.vector *= -1.0
		return temp

	def __mul__(self, other):

		assert isinstance(other, np.matrix) is True or \
		       isinstance(other, float)     is True or \
		       isinstance(other, complex)   is True or \
		       isinstance(other, int)       is True, "ERROR!  Type Error in LorentzVector rmul..."

		if isinstance(self, Momenta):
			if isinstance(other, np.matrix):
				temp = self.vector * other
				return Momenta(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=self.upper, incoming=self.incoming, physical=False)
			else:
				return Momenta(other * self.T(), other * self.X(), other * self.Y(), other * self.Z(), upper=self.upper, incoming=self.incoming, physical=False)
		else:
			if isinstance(other, np.matrix):
				temp = self.vector * other
				return LorentzVector(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=self.upper)
			else:
				return LorentzVector(other * self.T(), other * self.X(), other * self.Y(), other * self.Z(), upper=self.upper)

	def __pow__(self, exponent):

		assert isinstance(exponent, int) and exponent > 0, "ERROR! Wrong type of argument for exponent..." + str(type(exponent))

		# If even
		if exponent % 2 == 0:

			return self.dot(self) ** (exponent / 2)

		else:
			return (self.dot(self) ** ((exponent - 1) / 2)) * self

	def __getitem__(self, index):

		if index == 0:
			return self.T()
		if index == 1:
			return self.X()
		if index == 2:
			return self.Y()
		if index == 3:
			return self.Z()

		print "Wrong index for momenta access!"

	def __rmul__(self, other):

		assert isinstance(other, np.matrix) is True or \
		       isinstance(other, float)     is True or \
		       isinstance(other, complex)   is True or \
		       isinstance(other, int)       is True, "ERROR!  Type Error in LorentzVector rmul..."

		if isinstance(self, Momenta):
			if isinstance(other, np.matrix):
				temp = other * self.vector.transpose()
				return Momenta(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=self.upper, incoming=self.incoming, physical=False)
			else:
				return Momenta(other * self.T(), other * self.X(), other * self.Y(), other * self.Z(), upper=self.upper, incoming=self.incoming, physical=False)
		else:
			if isinstance(other, np.matrix):
				temp = other * self.vector.transpose()
				return LorentzVector(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=self.upper)
			else:
				return LorentzVector(other * self.T(), other * self.X(), other * self.Y(), other * self.Z(), upper=self.upper)

	def __eq__ (self, other):

		if other == None:
			return False

		if (self.T() - other.T()) > TOLERANCE:
			return False
		if (self.X() - other.X()) > TOLERANCE:
			return False
		if (self.Y() - other.Y()) > TOLERANCE:
			return False
		if (self.Z() - other.Z()) > TOLERANCE:
			return False

		return True

	# Dot product for two four vectors
	def dot(self, other):

		assert isinstance(other, LorentzVector) is True, "Warning! Dot product with non-FourVector type..." + str(type(other))

		# if self.index != other.index:
			# print "   -> Warning Lorentz indices not matching for contraction..."

		if other.upper == self.upper:
			s = -1.0
		else:
			s = +1.0

		return self.T() * other.T() + s * self.X() * other.X() + s * self.Y() * other.Y() + s * self.Z() * other.Z()

	def slashed(self):

		# Make sure the momenta has a lower index for computing slash
		if self.upper:
			temp = metric * self.vector.transpose()
		else:
			temp = self.vector.transpose()

		j = complex(0.0, 0.0)

		return temp.item(0) * gamma[0] + \
		       temp.item(1) * gamma[1] + \
		       temp.item(2) * gamma[2] + \
		       temp.item(3) * gamma[3]

	def __div__(self, other):

		if not (isinstance(other, float) or isinstance(other, complex) or isinstance(other, int)):
			print "ERROR!  Cant divide LorentzVector by non-float/complex/int type!"

		if isinstance(self, Momenta):
			return Momenta      (self.T() / other, self.X() / other, self.Y() / other, self.Z() / other, upper=self.upper, incoming=self.incoming, physical=False)
		else:
			return LorentzVector(self.T() / other, self.X() / other, self.Y() / other, self.Z() / other, self.upper)

	def Lower(self):

		# If we can lower the index
		if self.upper == True:

			temp = metric * self.vector

			return LorentzVector(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=False)

		else:
			print "Warning!  LorentzVector already has lower index..."

	def Raise(self):

		# If we cant raise the index
		if self.upper == True:
			print "Warning!  LorentzVector already has upper index..."

		else:
			temp = metric * self.vector

			return LorentzVector(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=True)

	# Is the LorentzVector covariant?
	def co(self):
		if self.upper:
			return False
		else:
			return True

	# Is the LorentzVector contravariant?
	def contra(self):
		if self.upper:
			return True
		else:
			return False

	def setVector(self, vector):
		self.vector = vector

	# Print the index for this current
	def printIndex(self):
		print self.index

	# Change the index for an instance
	def setIndex(self, index):
		self.index = index

	# Perform a Lorentz boost
	def boost(vector, bX, bY, bZ):

		b2 = bX ** 2 + bY ** 2 + bZ ** 2
		g  = (1 / sqrt(1 - b2)).real

		# Define the Lorentz matrix
		boostMatrix = np.matrix([[g, -g * bX, -g * bY, -g * bZ],
			                 [-g * bX, 1 + (g - 1) * bX * bX / b2,     (g - 1) * bX * bY / b2,     (g - 1) * bX * bZ / b2],
		                         [-g * bY,     (g - 1) * bX * bY / b2, 1 + (g - 1) * bY * bY / b2,     (g - 1) * bY * bZ / b2],
		                         [-g * bZ,     (g - 1) * bX * bZ / b2,     (g - 1) * bY * bZ / b2, 1 + (g - 1) * bZ * bZ / b2]])

		# Calculate the new boost
		temp = boostMatrix * vector

		# Return the new LorentzVector
		return LorentzVector(temp.e(), temp.pX(), temp.pY(), temp.pZ(), vector.upper)

# Derived class for particle momenta
class Momenta(LorentzVector):

	pCount = 0
	iCount = 0

	iList = []
	oList = []

	def __init__ (self, E, pX, pY, pZ, incoming=False, upper=True, physical=True):

		if self.pCount == None:
			self.pCount = 0

		# Save the incoming/outgoing state
		self.incoming = incoming

		if physical:
			self.ID = self.pCount
			Momenta.pCount += 1
			if incoming:
				Momenta.iCount += 1
				if Momenta.iCount > 2:
					print "Warning! More than two incoming momenta..."

		else:
			self.ID = -1

		LorentzVector.__init__(self, E, pX, pY, pZ, upper=upper)

		if incoming:
			Momenta.iList.append(self)
		else:
			Momenta.oList.append(self)

	@classmethod
	def Conserved(cls):

		temp  = Momenta.iList[0]
		temp += Momenta.iList[1]

		for p in Momenta.oList:
			temp -= p

		if temp.mass() < TOLERANCE:
			return True
		else:
			return False

	def isIncoming(self):
		if self.incoming:
			return True
		else:
			return False

	def mass(self):

		self.m = sqrt(self.dot(self))

		# If we square root has calculated a tiny imag. part b mistake ignore it
		if self.m.imag != 0.0 and self.m.imag < 1e-5:
			self.m = self.m.real

		if self.m.real < 0.0 or self.m.imag != 0:
			# print "Warning!  Negative or complex mass for momentas", self.m
			return self.m.real
		else:
			return self.m.real

	def mass2(self):

		self.m = sqrt(self.dot(self))

		# If we square root has calculated a tiny imag. part b mistake ignore it
		if self.m.imag != 0.0 and self.m.imag < 1e-5:
			self.m = self.m.real ** 2
		if self.m.real < 0.0 or self.m.imag != 0:
			return self.m.real ** 2
		else:
			return self.m.real ** 2

	def isMassless(self):

		if self.mass() == 0.0:
			return True
		else:
			return False

	def plus(self):
		return self.T() + self.Z()

	def minus(self):
		return self.T() - self.Z()

	def perp(self):
		return self.X() + i * self.Y()

	def boost(vector, bX, bY, bZ):

		b2 = bX ** 2 + bY ** 2 + bZ ** 2
		g  = (1 / sqrt(1 - b2)).real

		# Define the Lorentz matrix
		boostMatrix = np.matrix([[ g, -g * bX, -g * bY, -g * bZ],
			                 [-g * bX, 1 + (g - 1) * bX * bX / b2,     (g - 1) * bX * bY / b2,     (g - 1) * bX * bZ / b2],
		                         [-g * bY,     (g - 1) * bX * bY / b2, 1 + (g - 1) * bY * bY / b2,     (g - 1) * bY * bZ / b2],
		                         [-g * bZ,     (g - 1) * bX * bZ / b2,     (g - 1) * bY * bZ / b2, 1 + (g - 1) * bZ * bZ / b2]])

		# Calculate the new boost
		temp = boostMatrix * vector

		# Return the new LorentzVector
		return Momenta(temp.T(), temp.X(), temp.Y(), temp.Z(), vector.upper, physical=False)

	def decompose(self):

		mass = self.mass()

		if mass <= 0:

			print "ERROR!  No need to decompose massless momenta..."
			return

		# We can always boost the massive momenta to this one
		restFrameMomenta = Momenta(mass, 0.0,  0.0,  0.0, physical=False)

		# This can in turn be written as a linear combination of...
		basisVec1 = Momenta(mass,  mass, 0.0, 0.0, physical=False)
		basisVec2 = Momenta(mass, -mass, 0.0, 0.0, physical=False)

		# Calculate the boost needed to move to the comoving frame
		bX = - self.X() / self.T()
		bY = - self.Y() / self.T()
		bZ = - self.Z() / self.T()

		# Return the boosted vectors
		v1 = basisVec1.boost(bX, bY, bZ) / 2.0
		v2 = basisVec2.boost(bX, bY, bZ) / 2.0

		# Set them as not-oncoming
		v1.incoming = False
		v2.incoming = False

		return v1, v2

# A functor class whick holds an item (a current or a LorentzVector, ...) for each polarisation in the problem
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

	# Zero(ish) checker to fix screwed up numerics where things are basically zero but not quite
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

	# Define functor behaviour
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

		# Return the appropriate type of object (multiplication sholudnt change cast!)
		if isinstance(self, Current):
			return Current(self.spin1, self.index, self.spin2, upper=self.upper, Dict=temp)

		elif isinstance(self, polDictionary):
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

		# Return the appropriate type of object (division sholudnt change cast!)
		if isinstance(self, Current):
			return Current(self.spin1, self.index, self.spin2, upper=self.upper, Dict=temp)

		elif isinstance(self, polDictionary):
			return polDictionary(temp)

	def dot(self, other):

		temp = {}

		assert isinstance(other, LorentzVector) or \
		       isinstance(other, polDictionary), "ERROR! Wrong args. in dot..." + str(type((other)))

		# For each polarisation
		for key in self.polarisations:


			if isinstance(other, LorentzVector):# or isinstance(other, Current.current):

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

class Current(polDictionary):

	@classmethod
	def initFromDict(cls, instance, Dict):

		temp = deepcopy(instance)
		temp.polarisations = Dict
		return temp

	def __init__(self, spin1, index, spin2, upper=True, Dict=None):

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

	class current(LorentzVector):

		# Keep track of all the indices in the calculatation
		indices = {}

		def __init__ (self, spin1, index, spin2, upper=True, Vector=None):

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

class Spinor(object):

	def __init__ (self, momentum, barred=False, fermion=True, ID=None):

		# Spinors are defined using contravariant momenta so raise the index if passed covarant one
		if momentum.co():
			momentum = momentum.Raise()

		if ID == None:
			self.ID = momentum.ID
		else:
			self.ID = ID

		# Save the momentum for the spinor and calculate its mass
		self.momentum = momentum
		self.barred   = barred
		self.fermion  = fermion
		self.mass     = momentum.mass()
		self.pols     = {}

		# If the momentum is massless we can treat fermions and anti-fermions in the normal manner
		if self.mass == 0.0:

			if fermion:
				self.pols["+"] = self.spinor(momentum, "+", barred=barred, ID=self.ID)
				self.pols["-"] = self.spinor(momentum, "-", barred=barred, ID=self.ID)
			else:
				self.pols["+"] = self.spinor(momentum, "-", barred=barred, ID=self.ID)
				self.pols["-"] = self.spinor(momentum, "+", barred=barred, ID=self.ID)

		# Else use the decomposition trick
		else:

			# Decompose the massive momenta
			k_1, k_2 = momentum.decompose()

			# Double check the momenta are massless
			assert k_1.isMassless() and k_2.isMassless(), "Warning! Momentum decomposition broken..."

			# Define the massless spinors
			u_k1 = Spinor(k_1)
			u_k2 = Spinor(k_2)

			if fermion:
				sign = +1.0
			else:
				sign = -1.0

			if not barred:

				# Calculate the spinor vector
				temp1 = (u_k2 // u_k1) / momentum.mass() * u_k2("+") + sign * u_k1("-")
				temp2 = (u_k2 ** u_k1) / momentum.mass() * u_k2("-") + sign * u_k1("+")

				# Set the spinors
				self.pols["+"] = self.spinor(momentum, "+", barred=barred, vector=temp1, ID=self.ID)
				self.pols["-"] = self.spinor(momentum, "-", barred=barred, vector=temp2, ID=self.ID)

			if barred:

				# Calculate the spinor vector
				temp1 = (u_k1 ** u_k2) / momentum.mass() * u_k2.bar("+") + sign * u_k1.bar("-")
				temp2 = (u_k1 // u_k2) / momentum.mass() * u_k2.bar("-") + sign * u_k1.bar("+")

				# Set the spinors
				self.pols["+"] = self.spinor(momentum, "+", barred=barred, vector=temp1, ID=self.ID)
				self.pols["-"] = self.spinor(momentum, "-", barred=barred, vector=temp2, ID=self.ID)

	# Make the object a functor
	def __call__(self, polarisation):

		if polarisation == "+":
			return self.pols["+"]
		elif polarisation == "-":
			return self.pols["-"]
		else:
			print "Wrong polarisation!"

	def __str__(self):

		return "+: " + self("+").__str__() + "\n-: " + self("-").__str__()

	def __ne__ (self, other):

		if self("+").all() != other("+").all() or self("-").all() != other("-").all():
			return True
		else:
			return False

	# Overload the floor division operator to calculate spinor brackets like [21]
	def __floordiv__(self, other):

		assert isinstance(other, Spinor) is True, "ERROR!  Computing spinor bracket [i j] with non-Spinor type."
		assert self.isMassless()         is True, "ERROR!  [i j] only defined for massless spinors."

		if self.fermion:
			if other.fermion:
				return (self.bar('+') * other('-')).item(0)
			else:
				return (self.bar('+') * other("+")).item(0)
		else:
			if other.fermion:
				return (self.bar('-') * other("-")).item(0)
			else:
				return (self.bar('-') * other("+")).item(0)

	# Overload the power operator to calculate spinor brackets like <21>
	def __pow__(self, other):

		assert isinstance(other, Spinor) is True, "ERROR!  Computing spinor bracket <i j> with non-Spinor type."
		assert self.isMassless()         is True, "ERROR!  <i j> only defined for massless spinors."

		if self.fermion:
			if other.fermion:
				return (self("-").bar() * other("+")).item(0)
			else:
				return (self("-").bar() * other("-")).item(0)

		else:
			if other.fermion:
				return (self("+").bar() * other("+")).item(0)
			else:
				return (self("+").bar() * other("-")).item(0)

	def getMomentum(self):
		return self.momentum

	def isBarred(self):
		if self.barred:
			return True
		else:
			return False

	def isMassless(self):
		if self.mass == 0.0:
			return True
		else:
			return False

	# Clever barred calculation which makes the barred spinors look like functors too...but they arent really
	def bar(self, pol=None):

		temp = Spinor(self.momentum, barred = not self.barred, fermion=self.fermion, ID=self.ID)

		if pol == "+":
			return temp.pols["+"]
		elif pol == "-":
			return temp.pols["-"]
		else:
			return temp

	def diracEqn(self):

		tempBool = True

		# For massless spinors
		if self.isMassless():

			# Positive Helicity
			if np.linalg.norm(self.momentum.slashed() * self.pols["+"]) > 1e-5:
				print "Failed DE with +ve helicity\n", self.momentum.slashed() * self.pols["+"]
				tempBool = False

			# Negative Helicity
			if np.linalg.norm(self.momentum.slashed() * self.pols["-"]) > 1e-5:
				print "Failed DE with -ve helicity\n", self.momentum.slashed() * self.pols["-"]
				tempBool = False

		# For massive spinors
		else:

			if self.fermion == True:
				sign = +1.0
			else:
				sign = -1.0

			# Positive Helicity
			if np.linalg.norm((self.momentum.slashed() * self("+") - sign * self.mass * self("+"))) > 1e-5:
				tempBool = False

			# Negative Helicity
			if np.linalg.norm((self.momentum.slashed() * self("-") + sign * self.mass * self("-"))) > 1e-5:
				tempBool = False

		return tempBool

	class spinor:

		def __init__ (self, momentum, pol, barred=False, vector=None, ID=None):

			# Save the polarisation, momentum for the spinor and calculate its mass
			self.momentum = momentum
			self.barred   = barred
			self.pol      = pol
			self.ID       = ID

			# Calculate the spinor components if not constructing with a given vector
			if vector == None:

				# If the momentum is beamline assume its incoming
				if momentum.isIncoming():

					if pol == "+":
						# If the incoming parton is moving in the forward direction
						if momentum.Z() > 0.0:
							self.vector = np.matrix([sqrt(momentum.plus()), 0.0, 0.0, 0.0], complex).transpose()

						# Else its backwards moving
						else:
							self.vector = np.matrix([0.0, -sqrt(momentum.minus()), 0.0, 0.0], complex).transpose()

					else:

						# If the incoming parton is moving in the forward direction
						if momentum.Z() > 0.0:
							self.vector = np.matrix([0.0, 0.0, 0.0, -sqrt(momentum.plus())], complex).transpose()

						# Else its backwards moving
						else:
							self.vector = np.matrix([0.0, 0.0, -sqrt(momentum.minus()), 0.0], complex).transpose()

				# Else if the momenta is outgoing
				else:

					if pol == "+":
						self.vector = np.matrix([sqrt(momentum.plus()), sqrt(momentum.minus()) * momentum.perp() / abs(momentum.perp()), 0.0, 0.0], complex).transpose()
					else:
						self.vector = np.matrix([0.0, 0.0, sqrt(momentum.minus()) * momentum.perp().conjugate() / abs(momentum.perp()), -sqrt(momentum.plus())], complex).transpose()

				if self.barred:
					self.vector = self.vector.transpose().conjugate() * g0

			# Otherwise use the vector we have been passed
			else:
				self.vector = vector

		def __str__(self):
			return str(self.vector)

		def getPolarisation(self):
			return self.pol

		def getVector(self):
			return self.vector

		def bar(self):
			return Spinor.spinor(self.momentum, self.pol, barred = not self.barred, ID=self.ID)

		def __neg__(self):

			return Spinor.spinor(self.momentum, self.pol, barred=self.barred, vector= -1.0 * self.vector, ID=self.ID)

		def __eq__(self, other):

			assert isinstance(other, Spinor.spinor) or isinstance(other, np.matrix), "Wrong argument to equate with spinor..." + str(type(other))

			if isinstance(other, Spinor.spinor):
				for i in range(0, 4):
					if abs(self.vector.item(i) - other.vector.item(i)) > TOLERANCE:
						return False
			else:
				for i in range(0, 4):
					if abs(self.vector.item(i) - other.item(i)) > TOLERANCE:
						return False

			return True

		# Addition
		def __add__(self, other):

			if isinstance(other, Spinor.spinor):
				return self.vector + other.vector
			elif isinstance(other, np.matrix):
				return self.vector + other
			elif isinstance(other, complex):
				return self.vector + other
			else:
				print "ERROR!  No addition rule for type", type(other)

		# Subtraction
		def __sub__(self, other):

			if isinstance(other, Spinor.spinor):
				return self.vector - other.vector
			elif isinstance(other, np.matrix):
				return self.vector - other
			elif isinstance(other, complex):
				return self.vector - other
			else:
				print "ERROR!  No addition rule for type", type(other)

		# Right Addition - set array priority to high to get correct behaviour
		__array_priority__ = 100
		def __radd__(self, other):

			if isinstance(other, Spinor.spinor):
				return self.vector + other.vector
			elif isinstance(other, np.matrix):
				return self.vector + other
			elif isinstance(other, complex):
				return self.vector + other
			else:
				print "ERROR!  No addition rule for type", type(other)

		__array_priority__ = 100
		def __rsub__(self, other):

			return -1.0 * self.__radd__(-1.0 * other)

		# Left multiplication i.e. self * other
		def __mul__ (self, other):

			# Multiplication by spinor
			if isinstance(other, Spinor.spinor):
				return self.vector * other.vector

			# Multiplication by matrix
			elif isinstance(other, np.matrix):
				return self.vector * other

			# Multiplication by complex or float
			elif isinstance(other, float) or isinstance(other, complex):
				return self.vector * other

			else:
				print "ERROR! Unknown multiplication rule for type", type(other)

		def __div__(self, other):

			assert isinstance(other, float) is True or isinstance(other, complex) is True, "ERROR!  Dividing spinor by non-float type..."

			return self.vector / other

		# Right multiplication i.e. self * other
		def __rmul__ (self, other):

			# Multiplication by spinor
			if isinstance(other, Spinor.spinor):
				return other.vector * self.vector

			# Multiplication by matrix
			elif isinstance(other, np.matrix):
				return other * self.vector

			# Multiplication by complex or float
			elif isinstance(other, float) or isinstance(other, complex):
				return other * self.vector

			else:
				print "ERROR! Unknown multiplication rule for type", type(other)

		# Define non-equality for spinors
		def __ne__ (self, other):

			if self.pol.all() != other.pol.all():
				return True
			else:
				return False

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
