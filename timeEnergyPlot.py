import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


fileName = "/data/user/enpaudel/triggerStudy/energyTime.txt"

simEnergy = []
simTime = []
with open(fileName,"r") as fh:
	for line in fh:
		simEnergy.append(float(line.split()[0]))
		simTime.append(float(line.split()[1]))

print(simEnergy)
print(simTime)
simTime = np.asarray(simTime)/(60.0*60.0) #converts time to hour
simEnergy = (np.asarray(simEnergy)+9)

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
	# ax.legend(fontsize=26)
	ax.set_xlabel(r'log energy [eV]',fontsize=26)
	ax.set_ylabel(r"time [hrs.]", fontsize=26)
	plt.savefig("../plots/energyTimeCost.pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotEnergyTime(simEnergy,simTime)