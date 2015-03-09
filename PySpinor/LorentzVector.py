from Common import *

## Parent class for anything with a Lorentz index
#  Will require full rewrite inheriting from new Tensor class
#  and likely with LorentzIndex class to be written
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

	## Allows vector to be called with specific index name given by string argument
	def __call__(self, index):

		assert isinstance(index, str), "LorentzVector can only be called with index..."

		temp = deepcopy(self)
		temp.index = index

		return temp

	def __str__(self):
		return "(" + str(self.T()) + "; " + str(self.X()) + ", " + str(self.Y()) + ", " + str(self.Z()) + ")"

	## Subtraction for co- and contra-variant vectors (resp)
	def __sub__(self, other):

		assert self.upper == other.upper, "ERROR!  Can't subtract LorentzVectors with co and contra indices!"

		assert isinstance(other, LorentzVector), "ERROR! Can only subtract LorentzVector from LorentzVector"
		
		return LorentzVector(self.T() - other.T(), self.X() - other.X(), self.Y() - other.Y(), self.Z() - other.Z(), self.upper)

	## Addition for co- and contra-variant vectors (resp)
	def __add__(self, other):

		assert self.upper == other.upper, "ERROR!  Cant add LorentzVectors with co and contra indices!"

		assert isinstance(other, LorentzVector), "ERROR! Can only add LorentzVector to LorentzVector"

		return LorentzVector(self.T() + other.T(), self.X() + other.X(), self.Y() + other.Y(), self.Z() + other.Z(), self.upper)

	## Negation of vector elements
	def __neg__ (self):

		temp = deepcopy(self)
		temp.vector *= -1.0
		return temp

	## Defines vector-scalar and vector-matrix multiplication 
	#  (numpy matrix input required for latter)
	def __mul__(self, other):

		assert isinstance(other, np.matrix) is True or \
			   isinstance(other, float)     is True or \
			   isinstance(other, complex)   is True or \
			   isinstance(other, int)       is True, "ERROR!  Type Error in LorentzVector rmul..."

		if isinstance(other, np.matrix):
			temp = self.vector * other
			return LorentzVector(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=self.upper)
		else:
			return LorentzVector(other * self.T(), other * self.X(), other * self.Y(), other * self.Z(), upper=self.upper)

	## Allows for powers of vectors through the vector self inner product
	def __pow__(self, exponent):

		assert isinstance(exponent, int) and exponent > 0, "ERROR! Wrong type of argument for exponent..." + str(type(exponent))

		# If even
		if exponent % 2 == 0:

			return self.dot(self) ** (exponent / 2)

		else:
			return (self.dot(self) ** ((exponent - 1) / 2)) * self

	## Returns particular element of vector specified by numerical index
	def __getitem__(self, index):

		if index == 0:
			return self.T()
		if index == 1:
			return self.X()
		if index == 2:
			return self.Y()
		if index == 3:
			return self.Z()

		print "Wrong index for compenent access!"

	## Scalar-vector and matrix-vector multiplication
	def __rmul__(self, other):

		assert isinstance(other, np.matrix) is True or \
			   isinstance(other, float)     is True or \
			   isinstance(other, complex)   is True or \
			   isinstance(other, int)       is True, "ERROR!  Type Error in LorentzVector rmul..."


		if isinstance(other, np.matrix):
			temp = other * self.vector.transpose()
			return LorentzVector(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=self.upper)
		else:
			return LorentzVector(other * self.T(), other * self.X(), other * self.Y(), other * self.Z(), upper=self.upper)

	## Comparison of float elements of a vector up to some degree of tolerance
	#  TOLERANCE defined in Common.py
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

	## Dot product for two four vectors
	def dot(self, other):

		assert isinstance(other, LorentzVector) is True, "Warning! Dot product with non-FourVector type..." + str(type(other))

		# if self.index != other.index:
			# print "   -> Warning Lorentz indices not matching for contraction..."

		if other.upper == self.upper:
			s = -1.0
		else:
			s = +1.0

		return self.T() * other.T() + s * self.X() * other.X() + s * self.Y() * other.Y() + s * self.Z() * other.Z()

	## Feynman slash notation
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

	## Division by scalar
	def __div__(self, other):

		assert (isinstance(other, float) or isinstance(other, complex) or isinstance(other, int)), "ERROR!  Cant divide LorentzVector by non-float/complex/int type!"

		return LorentzVector(self.T() / other, self.X() / other, self.Y() / other, self.Z() / other, self.upper)

	## Lower index
	def Lower(self):

		# If we can lower the index
		if self.upper == True:

			temp = metric * self.vector

			return LorentzVector(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=False)

		else:
			print "Warning!  LorentzVector already has lower index..."

	## Raise index
	def Raise(self):

		# If we cant raise the index
		if self.upper == True:
			print "Warning!  LorentzVector already has upper index..."

		else:
			temp = metric * self.vector

			return LorentzVector(temp.item(0), temp.item(1), temp.item(2), temp.item(3), upper=True)

	## Bool: Is the LorentzVector covariant?
	def co(self):
		if self.upper:
			return False
		else:
			return True

	## Bool: Is the LorentzVector contravariant?
	def contra(self):
		if self.upper:
			return True
		else:
			return False

	## Sets elements given 1-d numpy array
	def setVector(self, vector):
		self.vector = vector

	## Print the index for this current
	def printIndex(self):
		print self.index

	## Change the index for an instance
	def setIndex(self, index):
		self.index = index

	## Perform a Lorentz boost
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