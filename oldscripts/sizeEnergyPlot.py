#!/bin/env python3

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fileDir = "/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/"
fileList = sorted(glob.glob(fileDir+"*Proc*"))
print(int(fileList[0].split("/")[-1][3:9]))

class sizeEnergy(object):
	"""docstring for sizeEnergy"""
	def __init__(self,fpath):
		super(sizeEnergy, self).__init__()
		self.size = os.path.getsize(str(fpath))
		self.ID = int(fpath.split("/")[-1][3:9])
		self.energy = (139 + (self.ID % 30))/10.0
seList = []
for ifiles in fileList:
	se = sizeEnergy(ifiles)
	seList.append(se)

def averageSize(energy):
	sizes = []
	for se in seList:
		if (se.energy - energy) < 0.5:
			sizes.append(se.size)
	sizeMean = np.mean(sizes)
	sizeSD = np.std(sizes)
	return sizeMean, sizeSD

energyList = sorted(set([se.energy for se in seList]))

sizeMeanList = []
sizeSDList = []
for ienergy in energyList:
	sizeMean,sizeSD = averageSize(ienergy)
	sizeMeanList.append(sizeMean/(1024*1024))
	sizeSDList.append(sizeSD/(1024*1024))


def plotEnergySize(energy,sizeMeanList,sizeSDList):
	"""plot energy size"""
	fig = plt.figure(figsize = (8,5))
	gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
	ax = fig.add_subplot(gs[0])
	ax.errorbar(energy,sizeMeanList,yerr=sizeSDList,fmt=".", ms=10,c="orangered",alpha=1,label="unzip")
	ax.set_title(r"proton showers",fontsize=26)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=26)
	ax.grid(True)
	# ax.set_yscale("log")
	# ax.legend(fontsize=20)
	ax.set_xlabel(r'log energy [eV]',fontsize=26)
	ax.set_ylabel(r"size [MB]", fontsize=26)
	plt.savefig("../plots/energySizeCost.pdf",transparent=False,bbox_inches='tight')
	plt.close()
plotEnergySize(energyList,sizeMeanList,sizeSDList)




