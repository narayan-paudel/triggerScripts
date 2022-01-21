#!/usr/bin/env python

import tables
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap

import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetFiltProcTrigCount.hdf5", help='Input files after running plotTrigHist.py.')
args = parser.parse_args()


def Rdisk(energy):
	"""calculates radius of disc in m overwhich shower core spread"""
	# energy = int(np.log10(energy)*10)/10.0
	ebin = np.log10(energy)
	# print("unit check",ebin,energy)
	return 800+(int(ebin)-5)*300 + (2*(int(ebin)-5)//3)*300 + ((int(ebin)-5)//3)*300

def showerArea(energy):
	'''Calculates area [km2] overwhich core of corsika shower was shifted for resimulation
	'''
	R = Rdisk(energy)
	return np.pi*R**2 * 10**(-6)

# f = tables.open_file("../simFiles/trigChkTest/DAT000005GenDetFiltProcTrigCount.hdf5")
# print(f.root.OfflineIceTopSLCTankPulsesTotalTankHit)
# print(f.root.)

# f = tables.open_file("/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetFiltProcTrigCount.hdf5")
# print(f.root)
# print(f.root.SLCTankCharge)
# print(f.root.SLCTankCharge.cols.value[:])

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

# fh_SLCTank = pd.read_hdf("/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetFiltProcTrigCount.hdf5",key="OfflineIceTopSLCTankPulsesTotalCharge")
# fh_HLCTank = pd.read_hdf("/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetFiltProcTrigCount.hdf5",key="OfflineIceTopHLCTankPulsesTotalCharge")
fh_primary = pd.read_hdf(str(args.input[0]),key="MCPrimary")
# print("columns",fh_primary.columns)
zenithList = fh_primary['zenith'].values
print("zenith list",zenithList)
print("zenith min max",min(zenithList),max(zenithList))
energyList = fh_primary['energy'].values
coreX = fh_primary['x'].values
coreY = fh_primary['y'].values
print(energyList)
# print("primary",fh_primary.columns)
# print("primary",fh_primary['zenith'])
SLCTankCharge = pd.read_hdf(str(args.input[0]),key="OfflineIceTopSLCTankPulsesTotalCharge")["value"].values
HLCTankCharge = pd.read_hdf(str(args.input[0]),key="OfflineIceTopHLCTankPulsesTotalCharge")["value"].values

SLCTankHit = pd.read_hdf(str(args.input[0]),key="OfflineIceTopSLCTankPulsesTotalHit")["value"].values
HLCTankHit = pd.read_hdf(str(args.input[0]),key="OfflineIceTopHLCTankPulsesTotalHit")["value"].values

SLCVEMHit = pd.read_hdf(str(args.input[0]),key="OfflineIceTopSLCVEMPulsesTotalHit")["value"].values
HLCVEMHit = pd.read_hdf(str(args.input[0]),key="OfflineIceTopHLCVEMPulsesTotalHit")["value"].values

ITGlobalTrig = pd.read_hdf(str(args.input[0]),key="ITGlobalTriggered")["value"].values
ITSMTTrig = pd.read_hdf(str(args.input[0]),key="ITSMTTriggered")["value"].values

# print(fh_HLCTankHit.columns)
# print(fh.head())
# print(fh.columns)
# print(fh["value"].values)

print("length compare",len(zenithList),len(HLCTankCharge))

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
class zenChargeHLCSLC(object):
	"""docstring for zenChargeHLCSLC"""
	def __init__(self,x,y,zenith,energy,SMTTrig,GlobalTrig,chargeHLC,chargeSLC,tankHitHLC,tankHitSLC,VEMHitHLC,VEMHitSLC):
		super(zenChargeHLCSLC, self).__init__()
		self.x = x
		self.y = y
		self.zenith = zenith
		self.energy = energy
		self.SMTTrig = SMTTrig
		self.GlobalTrig = GlobalTrig
		self.chargeHLC = chargeHLC
		self.chargeSLC = chargeSLC
		self.tankHitHLC = tankHitHLC
		self.tankHitSLC = tankHitSLC
		self.VEMHitHLC = VEMHitHLC
		self.VEMHitSLC = VEMHitSLC

zenChargeList = []
for x,y,zen,energy,smtTrig,globalTrig,hlc,slc,hlcHit,slcHit,vemHlcHit,vemSlcHit in zip(coreX,coreY,zenithList,energyList,ITSMTTrig,ITGlobalTrig,HLCTankCharge,SLCTankCharge,HLCTankHit,SLCTankHit,HLCVEMHit,SLCVEMHit):
	zc = zenChargeHLCSLC(x,y,zen,energy,smtTrig,globalTrig,hlc,slc,hlcHit,slcHit,vemHlcHit,vemSlcHit)
	zenChargeList.append(zc)

# slcLim = 100
# zenChargeList = np.asarray([zc for zc in zenChargeList if zc.chargeSLC < slcLim]) #to remove events with very high slc charge

def triggerEfficiency(n_trig,n_total):
	print("n_trig,n_total",n_trig,n_total)
	if n_total != 0:
		return (n_trig/n_total)
	else:
		return 0

def effectiveArea(n_trig,n_total,area):
	return triggerEfficiency(n_trig,n_total)*area

energyBins = 10**np.linspace(5, 8.0, 31) #[5*10**5,10**6,5*10**6,10**7,5*10**7,10**8]

def plotDiskRadius(energyBins):
	diskR = [Rdisk(energy) for energy in energyBins]
	print("disk radius",set(diskR))
	diskArea = [showerArea(energy) for energy in energyBins]
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	ax.plot(energyBins,diskR,".",ls='-',lw = 2.5,alpha=1)
	# ax.plot(energyBins,diskArea,".",ls='-',lw = 2.5,alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"energy [GeV]", fontsize=24)
	ax.set_ylabel(r"disk radius [m]", fontsize=24)
	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/diskRadius.pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotDiskRadius(energyBins)

def energyEfficiency(energyBins,zclist,ratioType):
	smtTrigEfficiencyList = []
	globalTrigEfficiencyList = []
	for nbin,binStart in enumerate(energyBins[:-1]):
		smtTrigInBin = [zc.SMTTrig for zc in zclist if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1]]
		n_total = sum([1 for zc in zclist if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1]])
		totalSMTTrig = sum(smtTrigInBin)
		globalTrigInBin = [zc.GlobalTrig for zc in zclist if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1]]
		totalGlobalTrig = sum(globalTrigInBin)
		areaDisk = showerArea(binStart)
		if ratioType == "trigEff":
			smtTrigEfficiency = triggerEfficiency(totalSMTTrig,n_total)
			globalTrigEfficiency = triggerEfficiency(totalGlobalTrig,n_total)
		elif ratioType == "effArea":
			smtTrigEfficiency = effectiveArea(totalSMTTrig,n_total,areaDisk)
			globalTrigEfficiency = effectiveArea(totalGlobalTrig,n_total,areaDisk)
		smtTrigEfficiencyList.append(smtTrigEfficiency)
		globalTrigEfficiencyList.append(globalTrigEfficiency)
	return smtTrigEfficiencyList,globalTrigEfficiencyList


def plotTrigEfficiency(zenChargeList):
	'''
	plots trigger efficiency in different zenith bins
	'''
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
		zclist = [zc for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
		smtTrigEfficiencyList, globalTrigEfficiencyList = energyEfficiency(energyBins,zclist,"trigEff")
		ax.plot(energyBins[:-1],smtTrigEfficiencyList,".",ls='-',lw = 2.5,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"energy [GeV]", fontsize=24)
	ax.set_ylabel(r"trigger efficiency", fontsize=24)
	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=12)
	plt.savefig(plotFolder+"/trigEfficiency.pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotTrigEfficiency(zenChargeList)

def plotEffectiveArea(zenChargeList):
	'''
	plots trigger efficiency in different zenith bins
	'''
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
		zclist = [zc for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
		smtTrigEfficiencyList, globalTrigEfficiencyList = energyEfficiency(energyBins,zclist,"effArea")
		ax.plot(energyBins[:-1],smtTrigEfficiencyList,".",ls='-',lw = 2.5,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"energy [GeV]", fontsize=24)
	ax.set_ylabel(r"effective area [$km^{2}$]", fontsize=24)
	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=12)
	plt.savefig(plotFolder+"/trigEffectiveArea.pdf",transparent=False,bbox_inches='tight')
	plt.close()
def getCircle(radius):
	theta = np.linspace(0,2*np.pi,100)
	return radius*np.cos(theta),radius*np.sin(theta)

def plotCoreScatter(zenChargeList,suffix_E,suffix_Z,trig):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	x = [zc.x for zc in zenChargeList]
	y = [zc.y for zc in zenChargeList]
	if int(suffix_E) < 6:
		xCirc,yCirc = getCircle(800)
		ax.plot(xCirc,yCirc,'-',c="blue",label="r = 800 m")
	elif int(suffix_E) < 7:
		xCirc,yCirc = getCircle(800)
		ax.plot(xCirc,yCirc,'-',c="blue",label="r = 800 m")
		xCirc,yCirc = getCircle(1100)
		ax.plot(xCirc,yCirc,'-',c="yellow",label="r = 1100 m")
	elif int(suffix_E) < 8:
		xCirc,yCirc = getCircle(800)
		ax.plot(xCirc,yCirc,'-',c="blue",label="r = 800 m")
		xCirc,yCirc = getCircle(1100)
		ax.plot(xCirc,yCirc,'-',c="yellow",label="r = 1100 m")
		xCirc,yCirc = getCircle(1700)
		ax.plot(xCirc,yCirc,'-',c="orangered",label="r = 1700 m")
	elif int(suffix_E) >= 8:
		xCirc,yCirc = getCircle(800)
		ax.plot(xCirc,yCirc,'-',c="blue",label="r = 800 m")
		xCirc,yCirc = getCircle(1100)
		ax.plot(xCirc,yCirc,'-',c="yellow",label="r = 1100 m")
		xCirc,yCirc = getCircle(1700)
		ax.plot(xCirc,yCirc,'-',c="orangered",label="r = 1700 m")
		xCirc,yCirc = getCircle(2600)
		ax.plot(xCirc,yCirc,'-',c="purple",label="r = 2600 m")
	ax.scatter(x,y,s=10)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"x [m]", fontsize=24)
	ax.set_ylabel(r"y [m]", fontsize=24)
	# ax.set_xlim(0,100)
	# ax.set_ylim(0,100)
	ax.grid(True,alpha=0.2)
	ax.set_aspect("equal")
	if suffix_E != "":
		ax.set_title(r"log(energy/GeV) = {0:.1f}, $\theta$ = {1:.1f}$^{{\circ}}$".format(float(suffix_E),np.arcsin(np.sqrt(float(suffix_Z)))*180/np.pi),fontsize=24)
	# print("saving file",plotFolder+"/scatterSLCHit"+str(suffix)+".png")
	plt.legend(fontsize=12)
	plt.savefig(plotFolder+"/scatterSLCHit"+str(suffix_E)+"sin2Z"+str(suffix_Z)+"trig"+str(trig)+".png",transparent=False,bbox_inches='tight')
	plt.close()


# def plotCoreScatterEnergy(zenChargeList,trig):
# 	for nbin,binStart in enumerate(energyBins[:-1]):
# 		if trig == True:
# 			zenChargeListInBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1] and abs(zc.SMTTrig - np.float64(1.0))<0.1]
# 		elif trig == False:
# 			zenChargeListInBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1]]
# 		if len(zenChargeListInBin) > 1:
# 			plotCoreScatter(zenChargeListInBin,np.log10(binStart),trig)


# plotCoreScatterEnergy(zenChargeList,trig=True)
# plotCoreScatterEnergy(zenChargeList,trig=False)

def plotCoreScatterEnergyZenith(zenChargeList,trig):
	for nbin_E,binStart_E in enumerate(energyBins[:-1]):
		if trig == True:
			zenChargeListInEBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin_E] and zc.energy < energyBins[nbin_E+1] and abs(zc.SMTTrig - np.float64(1.0))<0.1]
		elif trig == False:
			zenChargeListInEBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin_E] and zc.energy < energyBins[nbin_E+1]]
		for nbin_Z,binStart_Z in enumerate(sin2ZenBins[:-1]):
			zenChargeListInEZBin = [zc for zc in zenChargeListInEBin if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin_Z] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin_Z+1])]
			if len(zenChargeListInEZBin) > 1:
				plotCoreScatter(zenChargeListInEZBin,np.log10(binStart_E),binStart_Z, trig)
# plotCoreScatterEnergyZenith(zenChargeList,trig=True)
# plotCoreScatterEnergyZenith(zenChargeList,trig=False)


def plotHLCHitCoreDist(zenChargeList,suffix_E):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	HLCHits = [zc.tankHitHLC for zc in zenChargeList]
	coreDist = [np.sqrt(zc.x**2+zc.y**2) for zc in zenChargeList]
	ax.scatter(coreDist, HLCHits ,s=10)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"core distance [m]", fontsize=24)
	ax.set_ylabel(r"HLC tank hits", fontsize=24)
	# ax.set_xlim(0,100)
	# ax.set_ylim(0,100)
	ax.grid(True,alpha=0.2)
	# ax.set_aspect("equal")
	if suffix_E != "":
		ax.set_title(r"log(energy/GeV) = {0:.1f}".format(float(suffix_E)),fontsize=24)
	# print("saving file",plotFolder+"/scatterSLCHit"+str(suffix)+".png")
	# plt.legend(fontsize=12)
	plt.savefig(plotFolder+"/scatterHLCHitCoreDist"+str(suffix_E)+".png",transparent=False,bbox_inches='tight')
	plt.close()

def plotCoreScatterEnergyZenith(zenChargeList):
	for nbin_E,binStart_E in enumerate(energyBins[:-1]):
		zenChargeListInEBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin_E] and zc.energy < energyBins[nbin_E+1]]
		if len(zenChargeListInEBin) > 1:
			plotHLCHitCoreDist(zenChargeListInEBin,np.log10(binStart_E))
plotCoreScatterEnergyZenith(zenChargeList)



# plotEffectiveArea(zenChargeList)


# def plotHistChargeZenBinTriggeredFull(zenChargeList,LCType):
# 	'''
# 	plots histogram of charge for zenith bins.
# 	'''
# 	zenChargeTrigList = [zc for zc in zenChargeList if zc.SMTTrig == 1]
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	if LCType == "HLC":
# 		hlcChargeList = [zc.chargeHLC for zc in zenChargeList]
# 		chargeBins = np.linspace(min(hlcChargeList),max(hlcChargeList),6400)
# 		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
# 			chargeInBin = [zc.chargeHLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
# 			chargeTrigInBin = [zc.chargeHLC for zc in zenChargeTrigList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
# 			n,bins,patches = ax.hist(chargeInBin ,bins=chargeBins,histtype="step",lw=1.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
# 			ax.hist(chargeTrigInBin ,bins=chargeBins,histtype="step",ls="--",color=patches[0].get_facecolor(),lw=1.5,alpha=1)
# 			ax.set_xlim(-0.4,50)
# 	elif LCType == "SLC":
# 		slcChargeList = [zc.chargeSLC for zc in zenChargeList]
# 		chargeBins = np.linspace(min(slcChargeList),max(slcChargeList),80)
# 		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
# 			chargeInBin = [zc.chargeSLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
# 			chargeTrigInBin = [zc.chargeSLC for zc in zenChargeTrigList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
# 			n,bins,patches = ax.hist(chargeInBin ,bins=chargeBins,histtype="step",lw=1.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
# 			ax.hist(chargeTrigInBin ,bins=chargeBins,histtype="step",ls="--",color=patches[0].get_facecolor(),lw=1.5,alpha=1)
# 		ax.set_xlim(-0.4,50)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
# 	ax.set_xlabel(r"q [vem]", fontsize=24)
# 	ax.set_ylabel(r"$count$", fontsize=24)
# 	ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	# ax.set_yscale('log')
# 	ax.set_ylim(0,5100)
# 	ax.grid(True,alpha=0.2)
# 	l1=ax.legend(fontsize=20)
# 	point_dash = mlines.Line2D([], [], linestyle='--',lw=1.5,color='black', marker='',markersize=5, label=r"trig only")
# 	l2 = ax.legend(handles=[point_dash],loc="lower right",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
# 	ax.add_artist(l1)
# 	ax.add_artist(l2)
# 	plt.savefig(plotFolder+"/histChargeTriggZenBins"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
# 	plt.close()


def plotHistChargeZenBinTriggered(zenChargeList,LCType,zoom,ymax,suffix):
	'''
	plots histogram of charge for zenith bins.
	'''
	zenChargeTrigList = [zc for zc in zenChargeList if zc.SMTTrig == 1]
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	qbinSize = 1.0
	if LCType == "HLC":
		hlcChargeList = [zc.chargeHLC for zc in zenChargeList]
		chargeBins = list(np.arange(0,41,1.0))
		# chargeBins = list(np.arange(0,max(hlcChargeList),1.0))
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			chargeInBin = [zc.chargeHLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			chargeTrigInBin = [zc.chargeHLC for zc in zenChargeTrigList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			n,bins,patches = ax.hist(chargeInBin ,bins=chargeBins,histtype="step",lw=1.5,label=r"$\sin ^2 \theta$ = {:.1f}-{:.1f}".format(sin2ZenBins[nbin],sin2ZenBins[nbin+1]),alpha=1)
			ax.hist(chargeTrigInBin ,bins=chargeBins,histtype="step",ls="--",color=patches[0].get_facecolor(),lw=1.5,alpha=1)
	elif LCType == "SLC":
		slcChargeList = [zc.chargeSLC for zc in zenChargeList]
		chargeBins = list(np.arange(0,41,1.0))
		# chargeBins = list(np.arange(0,max(slcChargeList),1.0))
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			chargeInBin = [zc.chargeSLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			chargeTrigInBin = [zc.chargeSLC for zc in zenChargeTrigList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			n,bins,patches = ax.hist(chargeInBin ,bins=chargeBins,histtype="step",lw=1.5,label=r"$\sin ^2 \theta$ = {:.1f}-{:.1f}".format(sin2ZenBins[nbin],sin2ZenBins[nbin+1]),alpha=1)
			ax.hist(chargeTrigInBin ,bins=chargeBins,histtype="step",ls="--",color=patches[0].get_facecolor(),lw=1.5,alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"$Q_{{tot,{0}}}$ [vem]".format(LCType), fontsize=24)
	ax.set_ylabel(r"$count$", fontsize=24)
	if suffix != "":
		ax.set_title("log(E/GeV)="+str(suffix),fontsize=24)
	ax.set_yscale('log')
	# if zoom == True:
		# ax.set_xlim(-0.4,40)
		# ax.set_ylim(0,5430)
	# elif zoom == False:
	# 	ax.set_yscale('log')
	ax.grid(True,alpha=0.2)
	ax.set_ylim(None,ymax)
	l1=ax.legend(fontsize=12)
	point_dash = mlines.Line2D([], [], linestyle='--',lw=1.5,color='black', marker='',markersize=5, label=r"trig only")
	l2 = ax.legend(handles=[point_dash],loc="upper center",fontsize=12,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
	ax.add_artist(l1)
	ax.add_artist(l2)
	plt.savefig(plotFolder+"/histChargeTriggZenBins"+str(LCType)+"Z"+str(zoom)+"E"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotHistChargeZenBinTriggered(zenChargeList,"HLC",zoom=True,suffix="",ymax=10**4.6)
# plotHistChargeZenBinTriggered(zenChargeList,"SLC",zoom=True,suffix="",ymax=10**4.6)
# plotHistChargeZenBinTriggered(zenChargeList,"HLC",zoom=False,suffix="",ymax=10**4.2)
# plotHistChargeZenBinTriggered(zenChargeList,"SLC",zoom=False,suffix="",ymax=10**4.2)

def plotHistChargeZenEnergy(zenChargeList,LCType):
	for nbin,binStart in enumerate(energyBins[:-1]):
		zenChargeListInBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1]]
		if len(zenChargeListInBin) > 1:
			plotHistChargeZenBinTriggered(zenChargeListInBin,LCType,zoom=True,suffix=np.log10(binStart),ymax=10**3.4)

# plotHistChargeZenEnergy(zenChargeList,"HLC")
# plotHistChargeZenEnergy(zenChargeList,"SLC")



def plotHistChargeZenBin(zenChargeList,LCType):
	'''
	plots histogram of charge for zenith bins.
	'''
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	if LCType == "HLC":
		hlcChargeList = [zc.chargeHLC for zc in zenChargeList]
		chargeBins = np.linspace(min(hlcChargeList),max(hlcChargeList),40)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			chargeInBin = [zc.chargeHLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			ax.hist(chargeInBin ,bins=chargeBins,histtype="step",lw=2.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
	elif LCType == "SLC":
		slcChargeList = [zc.chargeSLC for zc in zenChargeList]
		chargeBins = np.linspace(min(slcChargeList),max(slcChargeList),40)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			chargeInBin = [zc.chargeSLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			ax.hist(chargeInBin ,bins=chargeBins,histtype="step",lw=2.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"q [vem]", fontsize=24)
	ax.set_ylabel(r"$count$", fontsize=24)
	ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	ax.set_yscale('log')
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histChargeZenBins"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotHistChargeZenBin(zenChargeList,"HLC")
# plotHistChargeZenBin(zenChargeList,"SLC")


def plot2DHistChargeZenBin(zenChargeList,LCType,suffix=""):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	# zenithBins = [0.0,0.2,0.4,0.6,0.8,1.0]
	zenithBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
	zenithList = [np.square(np.sin(zc.zenith)) for zc in zenChargeList]
	# zenithBins = np.linspace(min(zenithList),max(zenithList),4)
	if suffix != "":
		cmax = 10**4.6
	else:
		cmax = 10**3.4
	if LCType == "HLC":
		hlcChargeList = [zc.chargeHLC for zc in zenChargeList]
		hlcChargeBins = np.linspace(0,40,41)
		# hlcChargeBins = list(np.arange(0.0,max(hlcChargeList),100))
		hist = ax.hist2d(zenithList,hlcChargeList,bins=(zenithBins,hlcChargeBins),norm=mpl.colors.LogNorm(),cmin=None,cmax=cmax,alpha=1)
	elif LCType == "SLC":
		slcChargeList = [zc.chargeSLC for zc in zenChargeList]
		print("slc charge min max",min(slcChargeList),max(slcChargeList))
		# slcChargeBins = np.linspace(min(slcChargeList),max(slcChargeList),40)
		slcChargeBins = np.linspace(0,40,41)
		print("slcchargebins",slcChargeBins)
		# slcChargeBins = list(np.arange(0.0,max(slcChargeList),1))
		hist = ax.hist2d(zenithList,slcChargeList,bins=(zenithBins,slcChargeBins),norm=mpl.colors.LogNorm(),cmin=None,cmax=cmax,alpha=1)
	cbar = fig.colorbar(hist[3],ax=ax,pad=0.02)
	cbar.set_label(r"count",fontsize=24)
	cbar.ax.tick_params(labelsize=24)
	ax.tick_params(axis='both',which='both', direction='in',labelsize=24,pad=9.0)
	ax.set_ylabel(r"$Q_{{tot,{0}}}$ [vem]".format(LCType), fontsize=24)
	ax.set_xlabel(r"$\sin ^2 \theta$", fontsize=24)
	if suffix != "":
		ax.set_title("log(E/GeV)="+str(suffix),fontsize=24)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/2DhistChargeZenBins"+str(LCType)+"E"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plot2DHistChargeZenBin(zenChargeList,"HLC","")
# plot2DHistChargeZenBin(zenChargeList,"SLC","")


def plot2DHistChargeZenEnergy(zenChargeList,LCType):
	for nbin,binStart in enumerate(energyBins[:-1]):
		zenChargeListInBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1]]
		if len(zenChargeListInBin) > 1:
			plot2DHistChargeZenBin(zenChargeListInBin,LCType,np.log10(binStart))

# plot2DHistChargeZenEnergy(zenChargeList,"HLC")
# plot2DHistChargeZenEnergy(zenChargeList,"SLC")
		
###################################################################################################################
#################################3for hit counts###################################################################
def plotHistHitZenBinTriggered(zenChargeList,LCType,zoom,suffix):
	'''
	plots histogram of charge for zenith bins.
	'''
	zenChargeTrigList = [zc for zc in zenChargeList if zc.SMTTrig == 1]
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	qbinSize = 1
	if LCType == "HLC":
		hlcHitList = [zc.tankHitHLC for zc in zenChargeList]
		hitBins = np.linspace(0,40,41)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			hitInBin = [zc.tankHitHLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			hitTrigInBin = [zc.tankHitHLC for zc in zenChargeTrigList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			n,bins,patches = ax.hist(hitInBin ,bins=hitBins,histtype="step",lw=1.5,label=r"$\sin ^2 \theta$ = {:.1f}-{:.1f}".format(sin2ZenBins[nbin],sin2ZenBins[nbin+1]),alpha=1)
			ax.hist(hitTrigInBin ,bins=hitBins,histtype="step",ls="--",color=patches[0].get_facecolor(),lw=1.5,alpha=0.7)
	elif LCType == "SLC":
		slcHitList = [zc.tankHitSLC for zc in zenChargeList]
		hitBins = np.linspace(0,40,41)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			hitInBin = [zc.tankHitSLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			hitTrigInBin = [zc.tankHitSLC for zc in zenChargeTrigList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			n,bins,patches = ax.hist(hitInBin ,bins=hitBins,histtype="step",lw=1.5,label=r"$\sin ^2 \theta$ = {:.1f}-{:.1f}".format(sin2ZenBins[nbin],sin2ZenBins[nbin+1]),alpha=1)
			ax.hist(hitTrigInBin ,bins=hitBins,histtype="step",ls="--",color=patches[0].get_facecolor(),lw=1.5,alpha=0.7)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"{} tank hits".format(str(LCType)), fontsize=24)
	ax.set_ylabel(r"$count$", fontsize=24)
	if suffix != "":
		ax.set_title("log(E/GeV)="+str(suffix),fontsize=24)
	ax.set_yscale('log')
	# if zoom == True:
	# 	ax.set_xlim(-0.4,50)
	# 	ax.set_ylim(0,5430)
	# elif zoom == False:
	# 	ax.set_yscale('log')
	# ax.set_xlim(-0.4,40)
	ax.grid(True,alpha=0.2)
	l1=ax.legend(fontsize=12)
	point_dash = mlines.Line2D([], [], linestyle='--',lw=1.5,color='black', marker='',markersize=5, label=r"trig only")
	l2 = ax.legend(handles=[point_dash],loc="upper center",fontsize=12,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
	ax.add_artist(l1)
	ax.add_artist(l2)
	plt.savefig(plotFolder+"/histTankHitsTriggZenBins"+str(LCType)+"E"+str(suffix)+"Z"+str(zoom)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotHistHitZenBinTriggered(zenChargeList,"HLC",zoom=True,suffix="")
# plotHistHitZenBinTriggered(zenChargeList,"SLC",zoom=True,suffix="")
# plotHistHitZenBinTriggered(zenChargeList,"HLC",zoom=False,suffix="")
# plotHistHitZenBinTriggered(zenChargeList,"SLC",zoom=False,suffix="")

def plotHistHitZenEnergy(zenChargeList,LCType):
	for nbin,binStart in enumerate(energyBins[:-1]):
		zenChargeListInBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1]]
		if len(zenChargeListInBin) > 1:
			plotHistHitZenBinTriggered(zenChargeListInBin,LCType,zoom=True,suffix=np.log10(binStart))

# plotHistHitZenEnergy(zenChargeList,"HLC")
# plotHistHitZenEnergy(zenChargeList,"SLC")






def plotHistVEMHitZenBinTriggered(zenChargeList,LCType,zoom):
	'''
	plots histogram of charge for zenith bins.
	'''
	zenChargeTrigList = [zc for zc in zenChargeList if zc.SMTTrig == 1]
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	qbinSize = 1
	if LCType == "HLC":
		hlcHitList = [zc.VEMHitHLC for zc in zenChargeList]
		hitBins = np.linspace(0,40,41)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			hitInBin = [zc.VEMHitHLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			hitTrigInBin = [zc.VEMHitHLC for zc in zenChargeTrigList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			n,bins,patches = ax.hist(hitInBin ,bins=hitBins,histtype="step",lw=1.5,label=r"$\sin ^2 \theta$ = {:.1f}-{:.1f}".format(sin2ZenBins[nbin],sin2ZenBins[nbin+1]),alpha=1)
			ax.hist(hitTrigInBin ,bins=hitBins,histtype="step",ls="--",color=patches[0].get_facecolor(),lw=1.5,alpha=0.7)
	elif LCType == "SLC":
		slcHitList = [zc.VEMHitSLC for zc in zenChargeList]
		hitBins = np.linspace(0,40,41)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			hitInBin = [zc.VEMHitSLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			hitTrigInBin = [zc.VEMHitSLC for zc in zenChargeTrigList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			n,bins,patches = ax.hist(hitInBin ,bins=hitBins,histtype="step",lw=1.5,label=r"$\sin ^2 \theta$ = {:.1f}-{:.1f}".format(sin2ZenBins[nbin],sin2ZenBins[nbin+1]),alpha=1)
			ax.hist(hitTrigInBin ,bins=hitBins,histtype="step",ls="--",color=patches[0].get_facecolor(),lw=1.5,alpha=0.7)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"{} VEM hits".format(str(LCType)), fontsize=24)
	ax.set_ylabel(r"$count$", fontsize=24)
	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	ax.set_yscale('log')
	# if zoom == True:
	# 	ax.set_xlim(-0.4,50)
	# 	ax.set_ylim(0,5430)
	# elif zoom == False:
	# 	ax.set_yscale('log')
	# ax.set_xlim(-0.4,40)
	ax.grid(True,alpha=0.2)
	l1=ax.legend(fontsize=12)
	point_dash = mlines.Line2D([], [], linestyle='--',lw=1.5,color='black', marker='',markersize=5, label=r"trig only")
	l2 = ax.legend(handles=[point_dash],loc="upper center",fontsize=12,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
	ax.add_artist(l1)
	ax.add_artist(l2)
	plt.savefig(plotFolder+"/histVEMHitsTriggZenBins"+str(LCType)+"Z"+str(zoom)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotHistVEMHitZenBinTriggered(zenChargeList,"HLC",zoom=True)
# plotHistVEMHitZenBinTriggered(zenChargeList,"SLC",zoom=True)



def plotHistHitZenBin(zenChargeList,LCType):
	'''
	plots histogram of charge for zenith bins.
	'''
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	if LCType == "HLC":
		# hlcHitList = [zc.tankHitHLC for zc in zenChargeList]
		hlcHitList = HLCTankHit
		hitBins = np.linspace(0,40,41)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			hitInBin = [zc.tankHitHLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			ax.hist(hitInBin ,bins=hitBins,histtype="step",lw=2.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
	elif LCType == "SLC":
		# slcHitList = [zc.tankHitSLC for zc in zenChargeList]
		slcHitList = SLCTankHit
		hitBins = np.linspace(0,40,41)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			hitInBin = [zc.tankHitSLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			ax.hist(hitInBin ,bins=hitBins,histtype="step",lw=2.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"hits", fontsize=24)
	ax.set_ylabel(r"$count$", fontsize=24)
	ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	ax.set_yscale('log')
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histHitsZenBins"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()


def plotHist(value):
	'''
	plots histogram of charge for zenith bins.
	'''
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	hitBins = np.linspace(0,40,41)
	ax.hist(value ,bins=hitBins,histtype="step",lw=2.5,alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"hits", fontsize=24)
	ax.set_ylabel(r"$count$", fontsize=24)
	ax.set_yscale('log')
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histHitsSimple.pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotHist(HLCTankHit)
# print("HLC Hits",HLCTankHit)
# print("len hits",len(HLCTankHit))
# print("len charge",len(HLCTankCharge))
# print("HLC Hits")
# print("len hits",len(SLCTankHit))
# print("len charge",len(SLCTankCharge))



# def plot2DHistHitZenBin(zenChargeList,LCType):
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	# zenithBins = [0.0,0.2,0.4,0.6,0.8,1.0]
# 	zenithBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
# 	zenithList = [np.square(np.sin(zc.zenith)) for zc in zenChargeList]
# 	# zenithBins = np.linspace(min(zenithList),max(zenithList),4)
# 	if LCType == "HLC":
# 		hlcHitList = [zc.tankHitHLC for zc in zenChargeList]
# 		# hlcHitBins = list(np.arange(0.0,41.0,1.0))
# 		hlcHitBins = np.linspace(0,40,41)
# 		hist = ax.hist2d(zenithList,hlcHitList,bins=(50,50),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)#norm=mpl.colors.LogNorm(),,bins=(zenithBins,hlcHitBins)
# 	elif LCType == "SLC":
# 		slcHitList = [zc.tankHitSLC for zc in zenChargeList]
# 		# slcHitBins = list(np.arange(0.0,41.0,1.0))
# 		slcHitBins = np.linspace(0,40,41)
# 		hist = ax.hist2d(zenithList,slcHitList,bins=(50,50),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)#,bins=(zenithBins,slcHitBins)
# 	cbar = fig.colorbar(hist[3],ax=ax,pad=0.02)
# 	cbar.set_label(r"count",fontsize=24)
# 	cbar.ax.tick_params(labelsize=24)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=24,pad=9.0)
# 	ax.set_ylabel(r"{} hits".format(str(LCType)), fontsize=24)
# 	ax.set_xlabel(r"$\sin ^2 \theta$", fontsize=24)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	# ax.set_xscale('log')
# 	ax.grid(True,alpha=0.2)
# 	# ax.legend(fontsize=20)
# 	plt.savefig(plotFolder+"/2DhistHitZenBins"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# plot2DHistHitZenBin(zenChargeList,"HLC")
# plot2DHistHitZenBin(zenChargeList,"SLC")

def plot2DHistHitZenBin(zenChargeList,LCType,suffix=""):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))
	# zenithBins = [0.0,0.2,0.4,0.6,0.8,1.0]
	zenithBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
	zenithList = [np.square(np.sin(zc.zenith)) for zc in zenChargeList]
	# zenithBins = np.linspace(min(zenithList),max(zenithList),4)
	if LCType == "HLC":
		hlcHitList = [zc.tankHitHLC for zc in zenChargeList]
		# hlcHitBins = list(np.arange(0.0,41.0,1.0))
		hlcHitBins = np.linspace(0,120,121)
		H,xedges,yedges = np.histogram2d(zenithList,hlcHitList,bins=(zenithBins,hlcHitBins))
		H = H.T
		im = ax.imshow(H, interpolation="nearest",origin="lower",aspect="auto",norm=mpl.colors.LogNorm(),extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
		# hist = ax.hist2d(zenithList,hlcHitList,bins=(50,50),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)#norm=mpl.colors.LogNorm(),,bins=(zenithBins,hlcHitBins)
	elif LCType == "SLC":
		slcHitList = [zc.tankHitSLC for zc in zenChargeList]
		# slcHitBins = list(np.arange(0.0,41.0,1.0))
		slcHitBins = np.linspace(0,120,121)
		H,xedges,yedges = np.histogram2d(zenithList,slcHitList,bins=(zenithBins,slcHitBins))
		H = H.T
		im = ax.imshow(H, interpolation="nearest",origin="lower",aspect="auto",norm=mpl.colors.LogNorm(),extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
		# hist = ax.hist2d(zenithList,slcHitList,bins=(50,50),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)#,bins=(zenithBins,slcHitBins)
	cbar = fig.colorbar(im,ax=ax,pad=0.02)
	cbar.set_label(r"count",fontsize=24)
	cbar.ax.tick_params(labelsize=24)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24,pad=9.0)
	ax.set_ylabel(r"{} hits".format(str(LCType)), fontsize=24)
	ax.set_xlabel(r"$\sin ^2 \theta$", fontsize=24)
	if suffix != "":
		ax.set_title("log(E/GeV)="+str(suffix),fontsize=24)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/2DhistHitZenBins"+str(LCType)+"E"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plot2DHistHitZenBin(zenChargeList,"HLC","")
# plot2DHistHitZenBin(zenChargeList,"SLC","")

def plot2DHistHitZenEnergy(zenChargeList,LCType):
	for nbin,binStart in enumerate(energyBins[:-1]):
		zenChargeListInBin = [zc for zc in zenChargeList if zc.energy >= energyBins[nbin] and zc.energy < energyBins[nbin+1]]
		if len(zenChargeListInBin) > 1:
			plot2DHistHitZenBin(zenChargeListInBin,LCType,np.log10(binStart))

# plot2DHistHitZenEnergy(zenChargeList,"HLC")
# plot2DHistHitZenEnergy(zenChargeList,"SLC")







def plot2DHistHitZenTrig(zenChargeList,LCType):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))
	# zenithBins = [0.0,0.2,0.4,0.6,0.8,1.0]
	zenithBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
	zenithList = [np.square(np.sin(zc.zenith)) for zc in zenChargeList if abs(zc.SMTTrig - np.float64(1.0))<0.1]
	# zenithBins = np.linspace(min(zenithList),max(zenithList),4)
	if LCType == "HLC":
		hlcHitList = [zc.tankHitHLC for zc in zenChargeList if abs(zc.SMTTrig - np.float64(1.0))<0.1]
		# hlcHitBins = list(np.arange(0.0,41.0,1.0))
		hlcHitBins = np.linspace(0,40,41)
		H,xedges,yedges = np.histogram2d(zenithList,hlcHitList,bins=(zenithBins,hlcHitBins))
		H = H.T
		im = ax.imshow(H, interpolation="nearest",origin="lower",aspect="auto",norm=mpl.colors.LogNorm(),extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
		# hist = ax.hist2d(zenithList,hlcHitList,bins=(50,50),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)#norm=mpl.colors.LogNorm(),,bins=(zenithBins,hlcHitBins)
	elif LCType == "SLC":
		slcHitList = [zc.tankHitSLC for zc in zenChargeList if abs(zc.SMTTrig - np.float64(1.0))<0.1]
		# slcHitBins = list(np.arange(0.0,41.0,1.0))
		slcHitBins = np.linspace(0,35,36)
		H,xedges,yedges = np.histogram2d(zenithList,slcHitList,bins=(zenithBins,slcHitBins))
		H = H.T
		im = ax.imshow(H, interpolation="nearest",origin="lower",aspect="auto",norm=mpl.colors.LogNorm(),extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
		# hist = ax.hist2d(zenithList,slcHitList,bins=(50,50),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)#,bins=(zenithBins,slcHitBins)
	cbar = fig.colorbar(im,ax=ax,pad=0.02)
	cbar.set_label(r"count",fontsize=24)
	cbar.ax.tick_params(labelsize=24)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24,pad=9.0)
	ax.set_ylabel(r"{} hits".format(str(LCType)), fontsize=24)
	ax.set_xlabel(r"$\sin ^2 \theta$", fontsize=24)
	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/2DhistHitZenTrigs"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plot2DHistHitZenTrig(zenChargeList,"HLC")
# plot2DHistHitZenTrig(zenChargeList,"SLC")

# def plotHist(zenChargeList,LCType):
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))
# 	# zenithBins = [0.0,0.2,0.4,0.6,0.8,1.0]
# 	zenithBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
# 	zenithList = [np.square(np.sin(zc.zenith)) for zc in zenChargeList]
# 	# zenithBins = np.linspace(min(zenithList),max(zenithList),4)
# 	if LCType == "HLC":
# 		hlcHitList = [zc.tankHitHLC for zc in zenChargeList]
# 		# hlcHitBins = list(np.arange(0.0,41.0,1.0))
# 		hlcHitBins = np.linspace(0,40,41)
# 		ax.hist(hlcHitList)
# 		# hist = ax.hist2d(zenithList,hlcHitList,bins=(50,50),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)#norm=mpl.colors.LogNorm(),,bins=(zenithBins,hlcHitBins)
# 	elif LCType == "SLC":
# 		slcHitList = [zc.tankHitSLC for zc in zenChargeList]
# 		# slcHitBins = list(np.arange(0.0,41.0,1.0))
# 		slcHitBins = np.linspace(0,40,41)
# 		ax.hist(slcHitList)
# 		# hist = ax.hist2d(zenithList,slcHitList,bins=(50,50),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)#,bins=(zenithBins,slcHitBins)
# 	# cbar = fig.colorbar(hist[3],ax=ax,pad=0.02)
# 	# cbar.set_label(r"count",fontsize=24)
# 	# cbar.ax.tick_params(labelsize=24)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=24,pad=9.0)
# 	ax.set_ylabel(r"{} hits".format(str(LCType)), fontsize=24)
# 	ax.set_xlabel(r"$\sin ^2 \theta$", fontsize=24)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	# ax.set_xscale('log')
# 	ax.grid(True,alpha=0.2)
# 	# ax.legend(fontsize=20)
# 	plt.savefig(plotFolder+"/histCheck"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# plot2DHistHitZenBin(zenChargeList,"HLC")
# plot2DHistHitZenBin(zenChargeList,"SLC")




















def plotScatterHitZen(zenChargeList,LCType):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	zenithList = [np.square(np.sin(zc.zenith)) for zc in zenChargeList]
	if LCType == "HLC":
		hlcHitList = [zc.tankHitHLC for zc in zenChargeList]
		# hlcHitBins = list(np.arange(0.0,41.0,1.0))
		hlcHitBins = np.linspace(0,40,41)
		ax.scatter(zenithList,hlcHitList,s=10)
	elif LCType == "SLC":
		slcHitList = [zc.tankHitSLC for zc in zenChargeList]
		# slcHitBins = list(np.arange(0.0,41.0,1.0))
		slcHitBins = np.linspace(0,40,41)
		ax.scatter(zenithList,slcHitList,s=10)
	# cbar = fig.colorbar(hist[3],ax=ax,pad=0.02)
	# cbar.set_label(r"count",fontsize=24)
	# cbar.ax.tick_params(labelsize=24)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24,pad=9.0)
	ax.set_ylabel(r"{} hits".format(str(LCType)), fontsize=24)
	ax.set_xlabel(r"$\sin ^2 \theta$", fontsize=24)
	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/ScatterHitZen"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotScatterHitZen(zenChargeList,"HLC")
# plotScatterHitZen(zenChargeList,"SLC")

##############################################################################################################
########################plot scatter plots####################################################################
def plotScatterHLCSLCCharge(zenChargeList,suffix):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	hlcChargeList = [zc.chargeHLC for zc in zenChargeList]
	slcChargeList = [zc.chargeSLC for zc in zenChargeList]
	ax.scatter(hlcChargeList,slcChargeList,s=10)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"$Q_{tot,HLC}$", fontsize=24)
	ax.set_ylabel(r"$Q_{tot,SLC}$", fontsize=24)
	ax.set_xlim(0,100)
	ax.set_ylim(0,100)
	ax.grid(True,alpha=0.2)
	ax.set_aspect("equal")
	if suffix != "":
		ax.set_title(r"$\sin ^{{2}} \theta = {0} $".format(suffix),fontsize=24)
	plt.savefig(plotFolder+"/scatterHLCSLCChargeSin2Z"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotScatterHLCSLCCharge(zenChargeList,"")


def plotScatterHLCSLCChargeZenithBins(zenChargeList):
	for nbin,binStart in enumerate(sin2ZenBins[:-1]):
		zclist = [zc for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
		plotScatterHLCSLCCharge(zclist,str(binStart))

# plotScatterHLCSLCChargeZenithBins(zenChargeList)


def plotScatterHLCSLCHit(zenChargeList,suffix):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	hlcChargeList = [zc.tankHitHLC for zc in zenChargeList]
	slcChargeList = [zc.tankHitSLC for zc in zenChargeList]
	ax.scatter(hlcChargeList,slcChargeList,s=10)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"$Q_{tot,HLC}$", fontsize=24)
	ax.set_ylabel(r"$Q_{tot,SLC}$", fontsize=24)
	ax.set_xlim(0,100)
	ax.set_ylim(0,100)
	ax.grid(True,alpha=0.2)
	ax.set_aspect("equal")
	if suffix != "":
		ax.set_title(r"$\sin ^{{2}} \theta = {0} $".format(suffix),fontsize=24)
	plt.savefig(plotFolder+"/scatterHLCSLCHitSin2Z"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotScatterHLCSLCCharge(zenChargeList,"")


def plotScatterHLCSLCHitZenithBins(zenChargeList):
	for nbin,binStart in enumerate(sin2ZenBins[:-1]):
		zclist = [zc for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
		plotScatterHLCSLCHit(zclist,str(binStart))

plotScatterHLCSLCHitZenithBins(zenChargeList)

def plot2DHistHLCSLC(zenChargeList):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	# ax.set_xlim(0,50)
	# ax.set_ylim(0,50)
	hlcChargeList = [zc.chargeHLC for zc in zenChargeList]
	slcChargeList = [zc.chargeSLC for zc in zenChargeList]
	# hlcBins = list(np.arange(0,max(hlcChargeList),1.0))
	hlcBins = np.linspace(0,40,41)
	# slcBins = list(np.arange(0,max(slcChargeList),1.0))
	slcBins = np.linspace(0,40,41)
	hist = ax.hist2d(hlcChargeList,slcChargeList,bins=(hlcBins,slcBins),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)
	cbar = fig.colorbar(hist[3],ax=ax,pad=0.02)
	cbar.set_label(r"count",fontsize=24)
	cbar.ax.tick_params(labelsize=24)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"$Q_{tot,HLC}$", fontsize=24)
	ax.set_ylabel(r"$Q_{tot,SLC}$", fontsize=24)
	# ax.grid(True,alpha=0.2)
	ax.set_aspect("equal")
	plt.savefig(plotFolder+"/2DHLCSLCCharge2.pdf",transparent=False,bbox_inches='tight')
	plt.close()
# plot2DHistHLCSLC(zenChargeList)







def plotHistCharge(charges,LCType):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	chargeMean = np.mean(charges)
	chargeSD = np.std(charges)
	minLim = min(charges)
	maxLim = max(charges)
	n_bins = np.linspace(minLim,maxLim,200)
	# n_bins = np.linspace(0,90,9)
	ax.hist(charges,bins=n_bins,histtype="step",lw=2.5,edgecolor='#2ca02c',label=str(LCType)+" hist",alpha=1)
	ax.axvline(x=chargeMean, ymin=0, ymax=1,c='#2ca02c',ls="-.",label=r"$\bar{q}$")
	ax.axvspan(chargeMean - chargeSD,chargeMean + chargeSD,ymin=0, ymax=1,color='#2ca02c',ls="-",label=r"$\sigma$",alpha=0.1)
	# ax.hist(thetaList,bins=n_bins,histtype="step",lw=2.5,fc='orange',label="ratio: 12 - 75",alpha=0.75)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
	ax.tick_params(axis='both',which='major', length=8,width=1.2)
	ax.tick_params(axis='both',which='minor', length=4)
	ax.set_yscale('log')
	ax.grid(True)
	# ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
	# ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
	ax.set_xlabel(r"q [vem]", fontsize=28)
	ax.set_ylabel(r"$count$", fontsize=22)
	ax.set_ylim(1,None)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histCharge"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotHistCharge(SLCTankCharge,LCType="SLC")
# plotHistCharge(HLCTankCharge,LCType="HLC")





# import tables
# import pandas



# CORSIKA_ID = "DAT059871"
# outputDir = "/home/enpaudel/icecube/triggerStudy/simFiles/"


# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetTrigCount.hdf5", help='Input files after counting triggers')
# args = parser.parse_args()



# f = tables.open_file("/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetTrigCount.hdf5")
# print(f.root)
# print("trig count", f.root)
# # print("trig count", trigCount.root.ITSMTTriggered.cols.Event[:])
# # trigCount



def plotHistCharge2(charges):
	chargeName = [k for k,v in locals().iteritems() if v == charges][0]
	if "SLC" in chargeName:
		LCType="SLC"
	if "HLC" in chargeName:
		LCType="HLC"
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	chargeMean = np.mean(charges)
	chargeSD = np.std(charges)
	minLim = min(charges)
	maxLim = max(charges)
	n_bins = np.linspace(minLim,maxLim,20)
	# n_bins = np.linspace(0,90,9)
	ax.hist(charges,bins=n_bins,histtype="step",lw=2.5,edgecolor='#2ca02c',label="hist",alpha=1)
	ax.axvline(x=chargeMean, ymin=0, ymax=1,c='#2ca02c',ls="-.",label=r"$\bar{q}$")
	ax.axvspan(chargeMean - chargeSD,chargeMean + chargeSD,ymin=0, ymax=1,color='#2ca02c',ls="-",label=r"$\sigma$",alpha=0.1)
	# ax.hist(thetaList,bins=n_bins,histtype="step",lw=2.5,fc='orange',label="ratio: 12 - 75",alpha=0.75)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
	ax.tick_params(axis='both',which='major', length=8,width=1.2)
	ax.tick_params(axis='both',which='minor', length=4)
	# ax.set_yscale('log')
	ax.grid(True)
	# ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
	# ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
	ax.set_xlabel(r"q [vem]", fontsize=28)
	ax.set_ylabel(r"$count$", fontsize=22)
	ax.set_ylim(1,None)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histCharge"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()