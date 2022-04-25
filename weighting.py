from icecube.weighting.weighting import icetop_mc_weights
from icecube.weighting.fluxes import GaisserH4a_IT

import numpy as np
from enum import Enum


class ParticleType(Enum):
	Gamma       =    1
	PPlus       =   14
	He4Nucleus  =  402
	N14Nucleus  = 1407
	O16Nucleus  = 1608
	Al27Nucleus = 2713
	Fe56Nucleus = 5626
	NuE         =   66
	NuEBar      =   67
	NuMu        =   68
	NuMuBar     =   69
	NuTau       =  133
	NuTauBar    =  134
	MuMinus     =    6
	MuPlus      =    5

class PDGCode(Enum):
	Gamma       =         22
	PPlus       =       2212
	He4Nucleus  = 1000020040
	N14Nucleus  = 1000070140
	O16Nucleus  = 1000080160
	Al27Nucleus = 1000130270
	Fe56Nucleus = 1000260560
	NuE         =         12
	NuEBar      =        -12
	NuMu        =         14
	NuMuBar     =        -14
	NuTau       =         16
	NuTauBar    =        -16
	MuMinus     =         13
	MuPlus      =        -13


class GetWeight(object):
	"""docstring for GetWeight"""
	def __init__(self):
		super(GetWeight, self).__init__()
		# self.weight_file = "/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json"
		self.weight_file = "/home/enpaudel/icecube/triggerStudy/simFiles/dataset_info200Uneven.json"
		self.flux = GaisserH4a_IT()
		self.dsetP = 12360
		self.dsetHe = 12630
		self.dsetO = 12631
		self.dsetFe = 12362
		# self.dsetP = 10410
		# self.dsetHe = 11663
		# self.dsetO = 12605
		# self.dsetFe = 10889

	def getGenerator(self,dset):
		"""
		returns generator for given dset
		"""
		return icetop_mc_weights(dset,dataset_file=self.weight_file)

	def combinedGenerator(self,nfilesP,nfilesHe,nfilesO,nfilesFe):
		genP = self.getGenerator(self.dsetP)
		genHe = self.getGenerator(self.dsetHe)
		genO = self.getGenerator(self.dsetO)
		genFe = self.getGenerator(self.dsetFe)
		# print(nfilesP,genP, nfilesHe,genHe,nfilesO,genO, nfilesFe,genFe)
		self.generator = nfilesP*genP + nfilesHe*genHe + nfilesO*genO + nfilesFe*genFe

	def getWeight(self,nfilesP,nfilesHe,nfilesO,nfilesFe,zenith,energy,ptype):
		'''
		returns weights from energy and particle type and number of files of each type
		'''
		self.combinedGenerator(nfilesP,nfilesHe,nfilesO,nfilesFe)
		# print("check",energy,ptype,np.cos(zenith))
		# print("check",self.flux(energy,ptype))
		# print("check",self.generator(energy,ptype,np.cos(zenith)))
		return self.flux(energy,ptype)/self.generator(energy,ptype,np.cos(zenith))




