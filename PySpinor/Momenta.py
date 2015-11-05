from Common import *
from LorentzVector import LorentzVector

# Derived class for particle momenta
class Momenta(LorentzVector):

	pCount = 0
	iCount = 0

	iList = []
	oList = []

	def __init__ (self, E, pX, pY, pZ, label='', incoming=False, upper=True, physical=True):

		if self.pCount == None:
			self.pCount = 0

		# Save the incoming/outgoing state
		self.incoming = incoming

		# Save the label for printing purposes
		self.label = label

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

	@classmethod
	def byConservation(cls):

		temp  = Momenta.iList[0]
		temp += Momenta.iList[1]

		for p in Momenta.oList:
			temp -= p

		Momenta.oList.append(temp)

		return temp

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

	def __add__(self, other):

		assert self.upper == other.upper, "ERROR!  Can't add momenta with co and contra indices!"

		assert isinstance(other, Momenta), "ERROR! Can only add momenta to momenta"

		return Momenta(self.T() + other.T(), self.X() + other.X(), self.Y() + other.Y(), self.Z() + other.Z(), upper=self.upper, incoming=self.incoming, physical=False)

	def __sub__(self, other):

		assert self.upper == other.upper, "ERROR!  Can't subtract momenta with co and contra indices!"

		assert isinstance(other, Momenta), "ERROR! Can only subtract momenta from momenta"

		return Momenta(self.T() - other.T(), self.X() - other.X(), self.Y() - other.Y(), self.Z() - other.Z(), upper=self.upper, incoming=self.incoming, physical=False)

	def __mul__(self, other):

		assert isinstance(other, np.matrix) is True or \
		       isinstance(other, float)     is True or \
		       isinstance(other, complex)   is True or \
		       isinstance(other, int)       is True, "ERROR!  Type Error in momenta mul..."

		if isinstance(other, np.matrix):
			temp = self.vector * other
			return Momenta(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=self.upper, incoming=self.incoming, physical=False)
		else:
			return Momenta(other * self.T(), other * self.X(), other * self.Y(), other * self.Z(), upper=self.upper, incoming=self.incoming, physical=False)

	def __rmul__(self, other):

		assert isinstance(other, np.matrix) is True or \
		       isinstance(other, float)     is True or \
		       isinstance(other, complex)   is True or \
		       isinstance(other, int)       is True, "ERROR!  Type Error in momenta rmul..."

		if isinstance(other, np.matrix):
			temp = other * self.vector.transpose()
			return Momenta(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=self.upper, incoming=self.incoming, physical=False)
		else:
			return Momenta(other * self.T(), other * self.X(), other * self.Y(), other * self.Z(), upper=self.upper, incoming=self.incoming, physical=False)

	def __div__(self, other):

		assert (isinstance(other, float) or isinstance(other, complex) or isinstance(other, int)), "ERROR!  Cant divide Momenta by non-float/complex/int type!"

		return Momenta(self.T() / other, self.X() / other, self.Y() / other, self.Z() / other, upper=self.upper, incoming=self.incoming, physical=False)









