from Common import *

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