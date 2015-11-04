__array_priority__ = 100

from LorentzVector import LorentzVector
from Momenta import Momenta

# Functions to find the kinematic invariant s_pq
def s(p, q):

	assert isinstance(p, LorentzVector) is True and \
	       isinstance(q, LorentzVector) is True, "ERROR! s called with massive (or non-) LorentzVector..."

	return (p + q) ** 2

def abs2(complexVar):

	if isinstance(complexVar, complex):
		return (complexVar * complexVar.conjugate()).real
	else:
		return (complexVar * complexVar.conjugate()).real()

def avgFactor():

	return  2 ** - (Momenta.pCount - Momenta.iCount)

