from icecube.weighting.weighting import icetop_mc_weights
from icecube.weighting.fluxes import GaisserH4a_IT

import numpy as np


class GetWeight(object):
	"""docstring for GetWeight"""
	def __init__(self):
		super(GetWeight, self).__init__()
		self.weight_file = "/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json"
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
		self.generator = nfilesP*genP + nfilesHe*genHe + nfilesO*genO + nfilesFe*genFe

	def getWeight(self,nfilesP,nfilesHe,nfilesO,nfilesFe,zenith,energy,ptype):
		'''
		returns weights from energy and particle type and number of files of each type
		'''
		self.combinedGenerator(nfilesP,nfilesHe,nfilesO,nfilesFe)
		return self.flux(energy,ptype)/self.generator(energy,ptype,np.cos(zenith))




