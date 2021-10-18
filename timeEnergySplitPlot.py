import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


fileName = "/data/user/enpaudel/triggerStudy/energyTimeMulti.txt"
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']




class timeEnergy(object):
	"""docstring for timeEnergy"""
	def __init__(self, energy,tzip,tgen,tdet,tfilt,tproc):
		super(timeEnergy, self).__init__()
		self.energy = energy
		self.tzip = tzip
		self.tgen = tgen
		self.tdet = tdet
		self.tfilt = tfilt
		self.tproc = tproc

simEnergy = []
teList = []
with open(fileName,"r") as fh:
	for line in fh:
		te = timeEnergy(float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4]),float(line.split()[5]))
		teList.append(te)

energyList = sorted(set([te.energy for te in teList]))
print(energyList)

def averageTime(energy):
	tzip = []
	tgen = []
	tdet = []
	tfilt = []
	tproc = []
	for te in teList:
		if abs(te.energy - energy) < 0.5:
			tzip.append(te.tzip)
			tgen.append(te.tgen)
			tdet.append(te.tdet)
			tfilt.append(te.tfilt)
			tproc.append(te.tproc)
	tzip_avg = np.mean(tzip)
	tzip_sd = np.std(tzip)
	tgen_avg = np.mean(tgen)
	tgen_sd = np.std(tgen)
	tdet_avg = np.mean(tdet)
	tdet_sd = np.std(tdet)
	tfilt_avg = np.mean(tfilt)
	tfilt_sd = np.std(tfilt)
	tproc_avg = np.mean(tproc)
	tproc_sd = np.std(tproc)
	return tzip_avg,tzip_sd,tgen_avg,tgen_sd,tdet_avg,tdet_sd,tfilt_avg,tfilt_sd,tproc_avg,tproc_sd

tzip_avgList = []
tzip_sdList =  []
tgen_avgList = []
tgen_sdList =  []
tdet_avgList = []
tdet_sdList =  []
tfilt_avgList = []
tfilt_sdList =  []
tproc_avgList = []
tproc_sdList =  []
for ienergy in energyList:
	tzip_avg,tzip_sd,tgen_avg,tgen_sd,tdet_avg,tdet_sd,tfilt_avg,tfilt_sd,tproc_avg,tproc_sd=averageTime(ienergy)
	tzip_avgList.append(tzip_avg)
	tzip_sdList.append(tzip_sd)
	tgen_avgList.append(tgen_avg)
	tgen_sdList.append(tgen_sd)
	tdet_avgList.append(tdet_avg)
	tdet_sdList.append(tdet_sd)
	tfilt_avgList.append(tfilt_avg)
	tfilt_sdList.append(tfilt_sd)
	tproc_avgList.append(tproc_avg)
	tproc_sdList.append(tproc_sd)



tzip_avgList = np.asarray(tzip_avgList)/(3600)
tzip_sdList = np.asarray(tzip_sdList)/(3600)
tgen_avgList = np.asarray(tgen_avgList)/(3600)
tgen_sdList = np.asarray(tgen_sdList)/(3600)
tdet_avgList = np.asarray(tdet_avgList)/(3600)
tdet_sdList = np.asarray(tdet_sdList)/(3600)
tfilt_avgList = np.asarray(tfilt_avgList)/(3600)
tfilt_sdList = np.asarray(tfilt_sdList)/(3600)
tproc_avgList = np.asarray(tproc_avgList)/(3600)
tproc_sdList = np.asarray(tproc_sdList)/(3600)

energyList = (np.asarray(energyList)+9)




def plotEnergyTimeSplit(energy,tzip_avgList,tzip_sdList,tgen_avgList,tgen_sdList,tdet_avgList,tdet_sdList,tfilt_avgList,tfilt_sdList,tproc_avgList,tproc_sdList):
	"""plot energy time"""
	print("len",len(energy),len(tproc_avgList),len(tproc_sdList))
	fig = plt.figure(figsize = (8,5))
	gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
	ax = fig.add_subplot(gs[0])
	ax.errorbar(energy,tzip_avgList,yerr=tzip_sdList,fmt=".", ms=10,c="orangered",alpha=1,label="unzip")
	ax.errorbar(energy,tgen_avgList,yerr=tgen_sdList,fmt=".", ms=10,c='#1f77b4',alpha=1,label="gen")
	ax.errorbar(energy,tdet_avgList,yerr=tdet_sdList,fmt=".", ms=10,c='#2ca02c',alpha=1,label="det")
	ax.errorbar(energy,tfilt_avgList,yerr=tfilt_sdList,fmt=".", ms=10,c='#d62728',alpha=1,label="filt")
	ax.errorbar(energy,tproc_avgList,yerr=tproc_sdList,fmt=".", ms=10,c='#8c564b',alpha=1,label="proc")
	ax.set_title(r"proton showers",fontsize=26)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=26)
	ax.grid(True)
	# ax.set_yscale("log")
	ax.legend(fontsize=20)
	ax.set_xlabel(r'log energy [eV]',fontsize=26)
	ax.set_ylabel(r"time [hrs.]", fontsize=26)
	plt.savefig("../plots/energyTimeSplitCost.pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotEnergyTimeSplit(energyList,tzip_avgList,tzip_sdList,tgen_avgList,tgen_sdList,tdet_avgList,tdet_sdList,tfilt_avgList,tfilt_sdList,tproc_avgList,tproc_sdList)


def plotEnergyTimeSplitLog(energy,tzip_avgList,tzip_sdList,tgen_avgList,tgen_sdList,tdet_avgList,tdet_sdList,tfilt_avgList,tfilt_sdList,tproc_avgList,tproc_sdList):
	"""plot energy time"""
	print("len",len(energy),len(tproc_avgList),len(tproc_sdList))
	fig = plt.figure(figsize = (8,5))
	gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
	ax = fig.add_subplot(gs[0])
	ax.plot(energy,tzip_avgList,".", ms=10,c="orangered",alpha=1,label="unzip")
	ax.plot(energy,tgen_avgList,".", ms=10,c='#1f77b4',alpha=1,label="gen")
	ax.plot(energy,tdet_avgList,".", ms=10,c='#2ca02c',alpha=1,label="det")
	ax.plot(energy,tfilt_avgList,".", ms=10,c='#d62728',alpha=1,label="filt")
	ax.plot(energy,tproc_avgList,".", ms=10,c='#8c564b',alpha=1,label="proc")
	ax.set_title(r"proton showers",fontsize=26)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=26)
	ax.grid(True)
	ax.set_yscale("log")
	ax.legend(fontsize=20)
	ax.set_xlabel(r'log energy [eV]',fontsize=26)
	ax.set_ylabel(r"time [hrs.]", fontsize=26)
	plt.savefig("../plots/energyTimeSplitCostLog.pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotEnergyTimeSplitLog(energyList,tzip_avgList,tzip_sdList,tgen_avgList,tgen_sdList,tdet_avgList,tdet_sdList,tfilt_avgList,tfilt_sdList,tproc_avgList,tproc_sdList)


def plotEnergyTimeSplitDiff(energy,tzip_avgList,tzip_sdList,tgen_avgList,tgen_sdList,tdet_avgList,tdet_sdList,tfilt_avgList,tfilt_sdList,tproc_avgList,tproc_sdList):
	"""plot energy time"""
	print("len",len(energy),len(tproc_avgList),len(tproc_sdList))
	fig = plt.figure(figsize = (8,5))
	gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
	ax = fig.add_subplot(gs[0])
	ax.plot(energy,tzip_avgList,".", ms=10,c="orangered",alpha=1,label="unzip")
	ax.plot(energy,tgen_avgList-tzip_avgList,".", ms=10,c='#1f77b4',alpha=1,label="gen")
	ax.plot(energy,tdet_avgList-tgen_avgList,".", ms=10,c='#2ca02c',alpha=1,label="det")
	ax.plot(energy,tfilt_avgList-tdet_avgList,".", ms=10,c='#d62728',alpha=1,label="filt")
	ax.plot(energy,tproc_avgList-tfilt_avgList,".", ms=10,c='#8c564b',alpha=1,label="proc")
	ax.set_title(r"proton showers",fontsize=26)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=26)
	ax.grid(True)
	# ax.set_yscale("log")
	ax.legend(fontsize=20)
	ax.set_xlabel(r'log energy [eV]',fontsize=26)
	ax.set_ylabel(r"time [hrs.]", fontsize=26)
	plt.savefig("../plots/energyTimeSplitCostDiff.pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotEnergyTimeSplitDiff(energyList,tzip_avgList,tzip_sdList,tgen_avgList,tgen_sdList,tdet_avgList,tdet_sdList,tfilt_avgList,tfilt_sdList,tproc_avgList,tproc_sdList)





	 


# simTime = np.asarray(simTime)/(60.0*60.0) #converts time to hour
# simEnergy = (np.asarray(simEnergy)+9)

def plotEnergyTime(energy,time):
	"""plot energy time"""
	fig = plt.figure(figsize = (8,5))
	gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
	ax = fig.add_subplot(gs[0])
	ax.plot(energy,time,".", ms=10,c="orangered",alpha=1)	
	ax.set_title(r"proton showers",fontsize=26)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=26)
	ax.grid(True)
	ax.set_yscale("log")
	ax.legend(fontsize=26)
	ax.set_xlabel(r'log energy [eV]',fontsize=26)
	ax.set_ylabel(r"time [hrs.]", fontsize=26)
	plt.savefig("../plots/energyTimeCost.pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotEnergyTime(simEnergy,simTime)

print(averageTime(6.9))