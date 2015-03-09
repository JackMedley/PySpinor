import numpy as np
from Common import *

## General tensor class
#
# General tensor containing an n-dim numpy array
# for an n-rank tensor with a n-length list of indices
# which can be "upper" or "lower" (contra- or covariant, resp.)
class Tensor(object):

	## Tensor contructor
	def __init__(self, mat, iList, ulList):
		# super(Tensor, self).__init__()

		# Assert that matrix imput is numpy array type
		assert isinstance(mat, np.ndarray), "Error: array of elements must be numpy array type!"

		# Assert rank of tensor matches size of array
		assert 4**len(iList) == mat.size, "Error: rank of tensor does not match size of array!"

		# Assert that upper/lower list is same length as list of indices
		assert len(ulList) == len(iList), "Error: upper/lower list needs to be same length as index list!"

		self.tensorElements = mat
		self.indexList      = iList
		self.upperlowerList = ulList
		self.rank           = len(iList)

		# Generates a dictionary to return wether a given index is upper or lower
		self.upperlowerDict = {self.indexList[i] : self.upperlowerList[i] for i in range(0, self.rank)}


	# Overloading

	## At present, printing a tensor prints only its element array
	def __str__(self):
		return self.tensorElements.__str__()

	## Overloading addition for tensors of equal rank, indices and upper/lower lists
	def __add__(self, other):

		# Assert that rank and type are the same
		assert isinstance(other, Tensor), "Error: can only add tensor to other tensors!"
		assert self.rank == other.rank  , "Error: can only add tensors of equal rank!"

		# Assert that indices match
		for i in range(0, self.rank):
			assert self.indexList[i] == other.indexList[i] and self.upperlowerList == other.upperlowerList, "Error: to add, tensor indices must match!"

		return Tensor((self.tensorElements + other.tensorElements), self.indexList, self.upperlowerList)

	## Overloading subtraction for tensors of equal rank, indices and upper/lower lists
	def __sub__(self, other):

		# Assert that rank and type are the same
		assert isinstance(other, Tensor), "Error: can only add tensor to other tensors!"
		assert self.rank == other.rank  , "Error: can only add tensors of equal rank!"

		# Assert that indices match
		for i in range(0, self.rank):
			assert self.indexList[i] == other.indexList[i] and self.upperlowerList == other.upperlowerList, "Error: to add, tensor indices must match!"

		return Tensor((self.tensorElements - other.tensorElements), self.indexList, self.upperlowerList)

	## Definining negation of tensor elements
	def __neg__ (self):

		temp = deepcopy(self)
		temp.tensorElements *= -1.0
		return temp

	## Defining division by scalars
	def __div__(self, other):
		assert isinstance(other, float)    or \
		       isinstance(other, complex)  or \
		       isinstance(other, int)     , "Error:  Type Error in tensor __div__!"

		if isinstance(other, int): print "Careful when dividing by ints! Look for rounding errors."

		return Tensor(self.tensorElements/other, self.indexList, self.upperlowerList)

	## Lowers (in the Minkowski space sense) the index specified by the string indx
	def LowerIndex(self, indx):
		# Need to 
		assert isinstance(indx, str), "Error: provide index name as string to lower!"
		assert self.upperlowerDict[indx] == 'upper', "Error: can't lower a lower index!"

		# metricTensor = self.__class__(metric, [idx, 'placeholder'], ['lower' , 'lower'])
		count = 2
		indDict = {}
		for index in self.indexList:
			if index == indx:
				indDict[index] = 1
				self.upperlowerList[self.indexList.index(index)] = 'lower'
				self.upperlowerDict[index] = 'lower'
			elif self.indexList != indx:
				indDict[index] = count
				count += 1

		self.tensorElements = np.einsum(metric, [0, 1], self.tensorElements, [indDict[x] for x in self.indexList])

	## Raises (in the Minkowski space sense) the index specified by the string indx
	def RaiseIndex(self, indx):
		# Need to 
		assert isinstance(indx, str), "Error: provide index name as string to lower!"
		assert self.upperlowerDict[indx] == 'lower', "Error: can't lower a lower index!"

		count = 2
		indDict = {}
		for index in self.indexList:
			if index == indx:
				indDict[index] = 1
				self.upperlowerList[self.indexList.index(index)] = 'upper'
				self.upperlowerDict[index] = 'upper'
			elif self.indexList != indx:
				indDict[index] = count
				count += 1

		self.tensorElements = np.einsum(metric, [0, 1], self.tensorElements, [indDict[x] for x in self.indexList])



	## Overloading multiplication to take care of contractions/tensor products
	#  Also handles tensor-scalar multiplication (along with __rmul__ for scalar-tensor)
	def __mul__(self, other):
		assert isinstance(other, Tensor)   or \
		       isinstance(other, float)    or \
		       isinstance(other, complex)  or \
		       isinstance(other, int)     , "Error:  Type Error in tensor __mul__!"

		contractedIndices = []
		indDict = {}
		count = 0
		newIndexList = []
		newUpperLowerList = []

		if isinstance(other, Tensor):
			for index in self.indexList:
				indDict[index] = count
				count += 1
				if index in other.indexList:
					contractedIndices.append(index)
					if (self.upperlowerDict[index] == 'upper' and other.upperlowerDict[index] == 'upper'):
						other.LowerIndex(index)
					if (self.upperlowerDict[index] == 'lower' and other.upperlowerDict[index] == 'lower'):
						other.RaiseIndex(index)
				else:
					newIndexList.append(index)
					newUpperLowerList.append(self.upperlowerDict[index])

			for index in other.indexList:
				if index not in contractedIndices:
					indDict[index] = count
					count +=1
					newIndexList.append(index)
					newUpperLowerList.append(other.upperlowerDict[index])

			newTensorElements = np.einsum(self.tensorElements, [indDict[x] for x in self.indexList], \
				                     other.tensorElements, [indDict[x] for x in other.indexList] )

			if newIndexList == []: return newTensorElements

			return Tensor(newTensorElements, newIndexList, newUpperLowerList)
		else:
			return Tensor(other*self.tensorElements, self.indexList, self.upperlowerList)

	## Defines scalar-tensor multiplication
	def __rmul__(self, other):
		assert isinstance(other, float)    or \
		       isinstance(other, complex)  or \
		       isinstance(other, int)     , "Error:  Type Error in tensor __mul__!"

		return self.__mul__(other)



	# Getters and Setters

	## Returns a list of all tensor elements
	def GetElements(self):
		return self.tensorElements

	## Returns full list of indices
	def GetIndices(self):
		return self.indexList

	## Returns the list of contra/co- variant states of each index
	def GetUpperLowerList(self):
		return self.upperlowerList

	## Returns the dictionary of contra/co- variant states of each index. Index name is the key
	def GetUpperLowerDict(self):
		return self.upperlowerDict

	## Returns particular element of tensor
	#  specified by list of ints specifying position
	def GetElement(self, li):
		
		# Assert that list of indices is of corrent length
		assert len(li) == self.rank, "Error: list of indices does not match rank of tensor"

		# assert that each element in li is an integer and then arrange into numpy index form
		indices = []
		for el in li:
			assert isinstance(el, int), "Error: To return element, all indices must be integer"
			indices.append([el])

		# self.tensorElements[li] returns a single element list and so [0] is added to return correct type
		return self.tensorElements[li][0]

	## Set elements of tensor by supplying numpy array
	def SetElements(self, mat):
		self.tensorElements = mat

	## Set single element of tensor by specifying position to be set and value
	def SetElement(self, pos, val):
		
		# Assert that each element in li is an integer and then arrange into numpy index form
		indices = []
		for el in pos:
			assert isinstance(el, int), "Error: To return element, all indices must be integer"
			indices.append([el])

		self.tensorElements[indices] = val

	## Replaces tensor's index list with iList
	def SetIndices(self, iList):
		assert len(iList) == self.rank, "Error: length of index list needs to be rank!"
		self.indexList = [x for x in iList]

	## Replaces tensor's upper/lower list with upper/lower list
	def SetUpperLowerList(self, ulList):
		assert len(ulList) == self.rank, "Error: length of upper/lower list needs to be rank!"
		self.upperlowerList = [x for x in ulList]
		self.GetUpperLowerDict = {self.indexList[i]: ulList[i] for i in range(0, self.rank)}

	# Other methods

	## Defining a transpose (reversing order of all indices) which we
	#  may want to restrict to only rank 2 tensors
	def transpose(self):
		self.indexList = reversed(self.indexList)


# ToDo:
# Want to allow tensors of general dimension though this would require rewriting 
# Raise/Lower and prods, etc. to take into acount different metrics
# Deal with swapping/replacing indices (and therefore corresponding rows and columns) and the corresponding a/symm properties etc.





