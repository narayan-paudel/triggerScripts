#!/usr/bin/env python3

import os
import glob
import subprocess

import tables
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
from matplotlib.ticker import AutoMinorLocator,MultipleLocator
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap

from customColors import qualitative_colors

import numpy as np

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

from weighting import GetWeight, ParticleType, PDGCode
from inclinedTriggerTools import *

basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTest/"
hdf5NullListP = sorted(glob.glob(basePath+"p*Clean*.hdf5"))
hdf5NullListHe = sorted(glob.glob(basePath+"He*Clean*.hdf5"))
hdf5NullListO = sorted(glob.glob(basePath+"O*Clean*.hdf5"))
hdf5NullListFe = sorted(glob.glob(basePath+"Fe*Clean*.hdf5"))
# hdf5NullList = np.concatenate((hdf5NullListP,hdf5NullListHe,hdf5NullListO,hdf5NullListFe))
hdf5NullList = np.concatenate((hdf5NullListP,hdf5NullListHe,hdf5NullListO,hdf5NullListFe))
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"


# colorsList = ['#9467bd', '#e377c2','#1f77b4','#2ca02c','#bcbd22','#ff7f0e','#8c564b','#7f7f7f','#17becf','#d62728',
# 			'#4477AA', '#332288', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
# 			'#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']
colorsList = ['#1f77b4','#ff7f0e','#2ca02c','#8c564b','#9467bd', '#e377c2','#bcbd22','#7f7f7f','#17becf','#d62728',
			'#4477AA', '#332288','#2ca02c', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
			'#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']

colorsCustom = qualitative_colors(10)
colorsIter = iter(colorsCustom)

nfilesP = nCorFiles(hdf5NullListP)
nfilesHe = nCorFiles(hdf5NullListHe)
nfilesO = nCorFiles(hdf5NullListO)
nfilesFe = nCorFiles(hdf5NullListFe)
print("no of files",nfilesP,nfilesHe,nfilesO,nfilesFe)

trigWindow = 10**(-6) # in ns

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
# sin2ZenBins = [0.0,0.822]

H4aWeightList=[]
HLC6_5000List=[]

for hdfFile in hdf5NullList:
	# print(getValue_(hdfFile,key="H4aWeight"))
	H4aWeightList += list(getValue_(hdfFile,key="H4aWeight"))
	HLC6_5000List += list(getValue_(hdfFile,key="HLC6_5000"))
print("before weighting",sum(HLC6_5000List))
HLC6_5000WeightedList = []
for iweight,ihlc6 in zip(H4aWeightList,HLC6_5000List):
	HLC6_5000WeightedList.append(iweight*ihlc6)
print("after weighting",sum(HLC6_5000WeightedList))

evtList = extractEvents(hdf5NullList)
print("event list before",len(evtList))
evtList_contained = containedEvents(evtList,640)
# evtList_contained = containedEvents(evtList,410)
total_events = sum([1 for ievt in evtList])
print("total events",total_events)
# triggered_events = sum([ievt.ITSMTTriggered for ievt in evtList])
# print("trigger rate",total_events,triggered_events,triggered_events/total_events)


# # weights = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zenith,energy,ptype)
# # weithtsOne = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zenith[-1],energy[-1],ptype[-1])
# # print("weightCompare",weights[-1],weithtsOne,len(weights))
# # adjustedWeights = []



energyBins = 10**np.linspace(5, 8.0, 31)
# energyBins = 10**np.linspace(5, 8.0, 7)
energyBinslgE = np.linspace(5.0,8.9,4000)
energyBinCenter = [5.1,6.1,7.1,8.1]
print("energy bins",energyBins)

plotRadiusEnergy(energyBinslgE)

def plotCoreScatter_(x,y,suffix,title):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	# ax.scatter(x,y,s=10,alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"x [m]", fontsize=24)
	ax.set_ylabel(r"y [m]", fontsize=24)
	xCirc,yCirc = getCircle(800)
	ax.plot(xCirc,yCirc,'-',c="purple",lw=3.0,label="r = 800 m")
	xCirc,yCirc = getCircle(1100)
	ax.plot(xCirc,yCirc,'-',c="blue",lw=3.0,label="r = 1100 m")
	xCirc,yCirc = getCircle(1700)
	ax.plot(xCirc,yCirc,'-',c="orange",lw=3.0,label="r = 1700 m")
	# xCirc,yCirc = getCircle(2600)
	# ax.plot(xCirc,yCirc,'-',c="yellow",lw=3.0,label="r = 2600 m")
	ax.scatter(x,y,s=10,alpha=1)
	# ax.set_xlim(0,100)
	# ax.set_ylim(0,100)
	ax.grid(True,alpha=0.2)
	ax.set_title(title,fontsize=16)
	ax.set_aspect("equal")
	plt.legend(fontsize=12)
	plt.savefig(plotFolder+"/coreScatter"+str(suffix)+".png",transparent=False,bbox_inches='tight')
	plt.close()

# def plotCoreScatter(hdfFileList):
# 	x,y = getCore(hdfFileList)
# 	zenList,ptypeList,energyList = getZenithTypeEnergy(hdfFileList)
# 	plotCoreScatter_(x,y,"all_shower","all energy")
# 	energyBins = 10**(np.linspace(5.0,8.0,31))
# 	print("energyBins",energyBins)
# 	for n,nEnergy in enumerate(energyBins[:-1]):
# 		xInBin = []
# 		yInBin = []
# 		for ix,iy,ienergy in zip(x,y,energyList):
# 			if ienergy >= energyBins[n] and ienergy < energyBins[n+1]:
# 				xInBin.append(ix) 
# 				yInBin.append(iy)
# 		print("nBins",n,energyBins[n],energyBins[n+1])
# 		plotCoreScatter_(xInBin,yInBin,r"{:.1f}".format(np.log10(nEnergy)),r"lg(E[GeV]):{0:.1f}-{1:.1f}".format(np.log10(energyBins[n]),np.log10(energyBins[n+1])))

def plotCoreScatterEnergy(evtList,energyLow,energyHigh,filtKey):
	x = [ievt.coreX for ievt in evtList]
	y = [ievt.coreY for ievt in evtList]
	zenithBins = [np.arcsin(np.sqrt(i)) for i in np.linspace(0.0,1.0,11)][:-3]
	zenithBins.append(np.deg2rad(65))
	print("zenith bins",[np.sin(i)**2 for i in zenithBins])
	for n,nZenith in enumerate(zenithBins[:-1]):
		xInBin = []
		yInBin = []
		for ix,iy,izen,ienergy in zip(x,y,zenList,energyList):
			if ienergy >= 10**energyLow and ienergy < 10**energyHigh and izen >= zenithBins[n] and izen < zenithBins[n+1]:
				xInBin.append(ix) 
				yInBin.append(iy)
		plotCoreScatter_(xInBin,yInBin,r"energy_{:.1f}Zen{:.1f}_filt{}".format(energyLow,np.sin(nZenith)**2,str(filtKey)),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(np.rad2deg(zenithBins[n]),np.rad2deg(zenithBins[n+1])))

# plotCoreScatterEnergy(evtList,6.0,6.1,"None")
# plotCoreScatterEnergy([ievt for ievt in evtList if abs(ievt.ITSMTTriggered-1)<0.01] ,6.0,6.1,"sta3")
# plotCoreScatterEnergy(evtList,7.0,7.1,"None")
# plotCoreScatterEnergy([ievt for ievt in evtList if abs(ievt.ITSMTTriggered-1)<0.01],7.0,7.1,"sta3")
# plotCoreScatterEnergy(evtList,5.9,6.0,"None")
# plotCoreScatterEnergy([ievt for ievt in evtList if abs(ievt.ITSMTTriggered-1)<0.01],6.9,7.0,"sta3")


def plotCoreScatter_(x,y,suffix,title):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	# ax.scatter(x,y,s=10,alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"x [m]", fontsize=24)
	ax.set_ylabel(r"y [m]", fontsize=24)
	xCirc,yCirc = getCircle(800)
	ax.plot(xCirc,yCirc,'-',c="purple",lw=3.0,label="r = 800 m")
	xCirc,yCirc = getCircle(1100)
	ax.plot(xCirc,yCirc,'-',c="blue",lw=3.0,label="r = 1100 m")
	xCirc,yCirc = getCircle(1700)
	ax.plot(xCirc,yCirc,'-',c="orange",lw=3.0,label="r = 1700 m")
	# xCirc,yCirc = getCircle(2600)
	# ax.plot(xCirc,yCirc,'-',c="yellow",lw=3.0,label="r = 2600 m")
	ax.scatter(x,y,s=0.05,alpha=1)
	# ax.set_xlim(0,100)
	# ax.set_ylim(0,100)
	ax.grid(True,alpha=0.2)
	ax.set_title(title,fontsize=30)
	ax.set_aspect("equal")
	plt.legend(fontsize=12)
	plt.savefig(plotFolder+"/coreScatter"+str(suffix)+".png",transparent=False,bbox_inches='tight')
	plt.close()


def plotScatterCore(evtList,triggerType):
	'''
	plots core distance of triggered shower
	'''
	energyBins = 10**(np.linspace(5.0,8.0,4))
	evtList = selectTriggered(evtList,triggerType)
	for ebin, ebinStart in enumerate(energyBins[:-1]):
		fig = plt.figure(figsize=(8,5))
		gs = gridspec.GridSpec(nrows=1,ncols=1)
		ax = fig.add_subplot(gs[0])
		colorIter = iter(colorsList)
		lowEdge_E = energyBins[ebin]
		highEdge_E = energyBins[ebin+1]
		evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
		ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((energyBins[ebin])*10**9),np.log10((energyBins[ebin+1])*10**9)))
		for nbin, binStart in enumerate(sin2ZenBins[:-1]):
			# colorIter = iter(colorsCustom)
			lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
			highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
			evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
			distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
			x = [ievt.coreX for ievt in evtZenBin]
			y = [ievt.coreY for ievt in evtZenBin]
			if len(distanceList)>2:
				xbins = np.linspace(min(distanceList),max(distanceList),200)
			else:
				xbins = np.linspace(0,1000,200)
			ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
				r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),lw=2.5,alpha=1)
			plotCoreScatter_(x,y,r"energy_{:.1f}Zen{:.1f}Trig{}".format(np.log10(lowEdge_E),np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180/np.pi,triggerType),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180/np.pi))
		ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
		ax.set_xlabel(r"core distance [m]", fontsize=22)
		ax.set_ylabel(r"count", fontsize=22)
		# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
		# ax.set_xscale('log')
		# ax.set_ylim(0,1300)
		ax.set_ylim(0,600)
		ax.set_xlim(0,1710)
		# ax.yaxis.set_minor_locator(MultipleLocator(100))
		# ax.xaxis.set_minor_locator(MultipleLocator(0.1))
		ax.grid(True,alpha=0.6)
		ax.legend(fontsize=12,loc='upper left')
		plt.savefig(plotFolder+"/distanceScattTrig"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(energyBins[ebin])+9),transparent=False,bbox_inches='tight')
		plt.close()
# plotScatterCore(evtList,"None")
plotScatterCore(evtList,"sta3")
plotScatterCore(evtList,"sta1")

def triggerEfficiency(n_trig,n_total):
	# print("n_trig,n_total",n_trig,n_total)
	if n_total != 0:
		return (n_trig/n_total)
	else:
		return 0

def effectiveArea(n_trig,n_total,area):
	return triggerEfficiency(n_trig,n_total)*area


def plotTrigEfficiency(evtList,energyBins,weighting,triggerType,containment):
	'''
	plots trigger efficiency in different zenith bins
	'''
	if containment == True:
		evtList = containedEvents(evtList,640)
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	# colorIter = iter(colorsCustom)
	colorIter = iter(colorsList)
	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
		lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
		highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
		evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
		energyList = []
		efficiencyList = []
		for ebin, ebinStart in enumerate(energyBins[:-1]):
			lowEdge_E = energyBins[ebin]
			highEdge_E = energyBins[ebin+1]
			evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
			# totalEvts = len(evtEBin)
			weights = [ievt.H4aWeight for ievt in evtEBin]
			totalEvts = len(evtEBin)
			# sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
			if triggerType == "sta3":
				triggerList = [ievt.ITSMTTriggered for ievt in evtEBin]
			elif triggerType == "sta1":
				triggerList = [ievt.STA1Trigger for ievt in evtEBin]
			elif triggerType == "slc3":
				triggerList = [ievt.slc3Trig for ievt in evtEBin]
			elif triggerType == "slc4":
				triggerList = [ievt.slc4Trig for ievt in evtEBin]
			elif triggerType == "slc5":
				triggerList = [ievt.slc5Trig for ievt in evtEBin]
			trigEff = triggerEfficiency(sum(triggerList),totalEvts)
			efficiencyList.append(trigEff)
			energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))		
		ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
	ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
	ax.set_ylabel(r"trigger efficiency", fontsize=22)
	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	# ax.set_xscale('log')
	ax.set_ylim(0,1.01)
	ax.set_xlim(14,17)
	# ax.yaxis.set_minor_locator(MultipleLocator(100))
	ax.xaxis.set_minor_locator(MultipleLocator(0.1))
	ax.grid(True,alpha=0.5)
	ax.legend(fontsize=12)
	plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"Efficiency.pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="sta3",containment=False)
