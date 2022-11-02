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


plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

# colorsList = ['#9467bd', '#e377c2','#1f77b4','#2ca02c','#bcbd22','#ff7f0e','#8c564b','#7f7f7f','#17becf','#d62728',
# 			'#4477AA', '#332288', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
# 			'#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']
colorsList = ['#1f77b4','#ff7f0e','#2ca02c','#8c564b','#9467bd', '#e377c2','#bcbd22','#7f7f7f','#17becf','#d62728',
			'#4477AA', '#332288','#2ca02c', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
			'#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']

colorsCustom = qualitative_colors(10)
colorsIter = iter(colorsCustom)

def nCorFiles(hdfFileList):
	nfiles = 0
	for ihdf in hdfFileList:
		dataT = pd.read_hdf(ihdf,key="I3EventHeader")
		nfiles += len(set(dataT["Run"].values))
	return nfiles


trigWindow = 10**(-6) # in ns

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
# sin2ZenBins = [0.0,0.822]

class CREvent(object):
	"""docstring for CREvent"""
	def __init__(self,runID,eventID,CRType,energy,zenith,coreX,coreY):
		super(CREvent, self).__init__()
		self.runID = runID
		self.eventID = eventID
		self.CRType = CRType
		self.energy = energy
		self.zenith = zenith
		self.coreX = coreX
		self.coreY = coreY

	def addWeight(self,H4aWeight):
		self.H4aWeight = H4aWeight

	def addTrigger(self,ITSMTTriggered,STA1Trigger):
		self.ITSMTTriggered = ITSMTTriggered
		self.STA1Trigger = STA1Trigger

	def addFilter(self,ITSTA5_filter,ITSTA3_filter):
		self.ITSTA5_filter = ITSTA5_filter
		self.ITSTA3_filter = ITSTA3_filter
	def addDeltaTHLC(self,deltaTHLC):
		self.deltaTHLC = deltaTHLC
	def addDeltaTSLC(self,deltaT3SLC,deltaT4SLC,deltaT5SLC,deltaT6SLC,deltaT7SLC,deltaT8SLC):
		self.deltaT3SLC = deltaT3SLC
		self.deltaT4SLC = deltaT4SLC
		self.deltaT5SLC = deltaT5SLC
		self.deltaT6SLC = deltaT6SLC
		self.deltaT7SLC = deltaT7SLC
		self.deltaT8SLC = deltaT8SLC
	def addDeltaTSLCSC(self,deltaT3SLCSC,deltaT4SLCSC,deltaT5SLCSC,deltaT6SLCSC,deltaT7SLCSC,deltaT8SLCSC):
		self.deltaT3SLCSC = deltaT3SLCSC
		self.deltaT4SLCSC = deltaT4SLCSC
		self.deltaT5SLCSC = deltaT5SLCSC
		self.deltaT6SLCSC = deltaT6SLCSC
		self.deltaT7SLCSC = deltaT7SLCSC
		self.deltaT8SLCSC = deltaT8SLCSC
	def addSLCTrigger(self,slc3Trig,slc4Trig,slc5Trig,slc6Trig,slc7Trig,slc8Trig):
		self.slc3Trig = slc3Trig
		self.slc4Trig = slc4Trig
		self.slc5Trig = slc5Trig
		self.slc6Trig = slc6Trig
		self.slc7Trig = slc7Trig
		self.slc8Trig = slc8Trig
	def addDeltaTTank(self,deltaT3Tank,deltaT4Tank,deltaT5Tank,deltaT6Tank,deltaT7Tank,deltaT8Tank):
		self.deltaT3Tank = deltaT3Tank
		self.deltaT4Tank = deltaT4Tank
		self.deltaT5Tank = deltaT5Tank
		self.deltaT6Tank = deltaT6Tank
		self.deltaT7Tank = deltaT7Tank
		self.deltaT8Tank = deltaT8Tank
	def addDeltaTTankSC(self,deltaT3TankSC,deltaT4TankSC,deltaT5TankSC,deltaT6TankSC,deltaT7TankSC,deltaT8TankSC):
		self.deltaT3TankSC = deltaT3TankSC
		self.deltaT4TankSC = deltaT4TankSC
		self.deltaT5TankSC = deltaT5TankSC
		self.deltaT6TankSC = deltaT6TankSC
		self.deltaT7TankSC = deltaT7TankSC
		self.deltaT8TankSC = deltaT8TankSC
	def addTankTrigger(self,tank1Trig,tank3Trig,tank4Trig,tank5Trig,tank6Trig,tank7Trig,tank8Trig):
		self.tank1Trig = tank1Trig
		self.tank3Trig = tank3Trig
		self.tank4Trig = tank4Trig
		self.tank5Trig = tank5Trig
		self.tank6Trig = tank6Trig
		self.tank7Trig = tank7Trig
		self.tank8Trig = tank8Trig
	def addSLCHLC(self,nSLC,nHLC):
		self.nSLC = nSLC
		self.nHLC = nHLC



def getValue_(hdfFile,key):
	'''extract value of given key in hdf file 
	as single array
	'''
	dataT = pd.read_hdf(hdfFile,key=key)
	return dataT["value"].values

def getValue(hdfFileList,key):
	hitList = np.array([])
	for ihdf in hdfFileList:
		ihit = getValue_(ihdf,key)
		hitList = np.concatenate((hitList,ihit))
	return hitList


def extractEvents_(hdfFile):
	"""
	extracts events from the hdfFile
	"""
	mcPrim_df = pd.read_hdf(hdfFile,key="MCPrimary")
	evtHead_df = pd.read_hdf(hdfFile,key="I3EventHeader")
	eventList = evtHead_df["Event"].values
	runList = evtHead_df["Run"].values
	primaryList = mcPrim_df["type"].values
	energyList = mcPrim_df["energy"].values
	zenithList = mcPrim_df["zenith"].values
	coreXList = mcPrim_df["x"].values
	coreYList = mcPrim_df["y"].values
	H4aWeightList = getValue_(hdfFile,key="H4aWeight")
	ITSMTTriggeredList = getValue_(hdfFile,key="ITSMTTriggered")
	STA1TriggerList = getValue_(hdfFile,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_isSTA1")
	ITSTA5_filterList = getValue_(hdfFile,key="IceTopSTA5_13_filter")
	ITSTA3_filterList = getValue_(hdfFile,key="SDST_IceTopSTA3_13_filter")
	hlcHitList = getValue_(hdfFile,key="OfflineIceTopHLCTankPulsesTotalHit")
	deltaTSLC3List = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_3SLCTankHit_t")
	deltaTSLC3ListSC = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_3SLCTankHit_tSC")
	deltaTSLC4List = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_4SLCTankHit_t")
	deltaTSLC4ListSC = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_4SLCTankHit_tSC")
	deltaTSLC5List = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_5SLCTankHit_t")
	deltaTSLC5ListSC = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_5SLCTankHit_tSC")
	deltaTSLC6List = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_6SLCTankHit_t")
	deltaTSLC6ListSC = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_6SLCTankHit_tSC")
	deltaTSLC7List = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_7SLCTankHit_t")
	deltaTSLC7ListSC = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_7SLCTankHit_tSC")
	deltaTSLC8List = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_8SLCTankHit_t")
	deltaTSLC8ListSC = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_8SLCTankHit_tSC")	
	slc3TrigList = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank3")
	slc4TrigList = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank4")
	slc5TrigList = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank5")
	slc6TrigList = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank6")
	slc7TrigList = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank7")
	slc8TrigList = getValue_(hdfFile,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank8")
	nSLC = getValue_(hdfFile,key="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalHit")
	nHLC = getValue_(hdfFile,key="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalHit")
	nStations = getValue_(hdfFile,key="HLCSLC_hitStations")
	nTanks = getValue_(hdfFile,key="HLCSLC_hitTanks")
	deltaT3List = getValue_(hdfFile,key="HLCSLC_3TankHit_t")
	deltaT4List = getValue_(hdfFile,key="HLCSLC_4TankHit_t")
	deltaT5List = getValue_(hdfFile,key="HLCSLC_5TankHit_t")
	deltaT6List = getValue_(hdfFile,key="HLCSLC_6TankHit_t")
	deltaT7List = getValue_(hdfFile,key="HLCSLC_7TankHit_t")
	deltaT8List = getValue_(hdfFile,key="HLCSLC_8TankHit_t")
	deltaT3ListSC = getValue_(hdfFile,key="HLCSLC_3TankHit_tSC")
	deltaT4ListSC = getValue_(hdfFile,key="HLCSLC_4TankHit_tSC")
	deltaT5ListSC = getValue_(hdfFile,key="HLCSLC_5TankHit_tSC")
	deltaT6ListSC = getValue_(hdfFile,key="HLCSLC_6TankHit_tSC")
	deltaT7ListSC = getValue_(hdfFile,key="HLCSLC_7TankHit_tSC")
	deltaT8ListSC = getValue_(hdfFile,key="HLCSLC_8TankHit_tSC")
	TankHit1 = getValue_(hdfFile,key="HLCSLC_isTank1")
	TankHit3 = getValue_(hdfFile,key="HLCSLC_isTank3")
	TankHit4 = getValue_(hdfFile,key="HLCSLC_isTank4")
	TankHit5 = getValue_(hdfFile,key="HLCSLC_isTank5")
	TankHit6 = getValue_(hdfFile,key="HLCSLC_isTank6")
	TankHit7 = getValue_(hdfFile,key="HLCSLC_isTank7")
	TankHit8 = getValue_(hdfFile,key="HLCSLC_isTank8")
	# deltaTHLCList = pd.read_hdf(hdfFile,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t")["item"].values
	evtObjList = []
	for nEvt,evtID in enumerate(eventList):
		# if STA1TriggerList[nEvt] > 0:
			# print("primary,runID,eventID,STA1,STA3,STA5",primaryList[nEvt],runList[nEvt],evtID,ITSMTTriggeredList[nEvt],STA1TriggerList[nEvt],ITSTA3_filterList[nEvt],ITSTA5_filterList[nEvt])
		thisEvt = CREvent(runList[nEvt],evtID,primaryList[nEvt],energyList[nEvt],zenithList[nEvt],coreXList[nEvt],coreYList[nEvt])
		thisEvt.addWeight(H4aWeightList[nEvt])
		thisEvt.addTrigger(ITSMTTriggeredList[nEvt],STA1TriggerList[nEvt])
		thisEvt.addFilter(ITSTA5_filterList[nEvt],ITSTA3_filterList[nEvt])
		thisEvt.hlcHit = hlcHitList[nEvt]
		thisEvt.addDeltaTSLC(deltaTSLC3List[nEvt],deltaTSLC4List[nEvt],deltaTSLC5List[nEvt],deltaTSLC6List[nEvt],deltaTSLC7List[nEvt],deltaTSLC8List[nEvt])
		thisEvt.addDeltaTSLCSC(deltaTSLC3ListSC[nEvt],deltaTSLC4ListSC[nEvt],deltaTSLC5ListSC[nEvt],deltaTSLC6ListSC[nEvt],deltaTSLC7ListSC[nEvt],deltaTSLC8ListSC[nEvt])
		thisEvt.addDeltaTTank(deltaT3List[nEvt],deltaT4List[nEvt],deltaT5List[nEvt],deltaT6List[nEvt],deltaT7List[nEvt],deltaT8List[nEvt])
		thisEvt.addDeltaTTankSC(deltaT3ListSC[nEvt],deltaT4ListSC[nEvt],deltaT5ListSC[nEvt],deltaT6ListSC[nEvt],deltaT7ListSC[nEvt],deltaT8ListSC[nEvt])
		thisEvt.addSLCTrigger(slc3TrigList[nEvt],slc4TrigList[nEvt],slc5TrigList[nEvt],slc6TrigList[nEvt],slc7TrigList[nEvt],slc8TrigList[nEvt])
		thisEvt.addTankTrigger(TankHit1[nEvt],TankHit3[nEvt],TankHit4[nEvt],TankHit5[nEvt],TankHit6[nEvt],TankHit7[nEvt],TankHit8[nEvt])
		thisEvt.addSLCHLC(nSLC[nEvt],nHLC[nEvt])
		thisEvt.nStations = nStations[nEvt]
		thisEvt.nTanks = nTanks[nEvt]
		# thisEvt.addDeltaTHLC(deltaTHLCList[nEvt])
		evtObjList.append(thisEvt)
	return evtObjList


def extractEvents(hdf5List):
	evtList = []
	for ihdf in hdf5List:
		evtList += extractEvents_(ihdf)
	return evtList

def removeFirstCore(evtList):
	return [ievt for ievt in evtList if abs(ievt.eventID % 100) > 0]

def extractEventsOfficial_(hdfFile):
	"""
	extracts events from the hdfFile
	"""
	mcPrim_df = pd.read_hdf(hdfFile,key="MCPrimary")
	evtHead_df = pd.read_hdf(hdfFile,key="I3EventHeader")
	eventList = evtHead_df["Event"].values
	runList = evtHead_df["Run"].values
	primaryList = mcPrim_df["type"].values
	energyList = mcPrim_df["energy"].values
	zenithList = mcPrim_df["zenith"].values
	coreXList = mcPrim_df["x"].values
	coreYList = mcPrim_df["y"].values
	H4aWeightList = getValue_(hdfFile,key="H4aWeight")
	# ITSMTTriggeredList = getValue_(hdfFile,key="ITSMTTriggered")
	# STA1TriggerList = getValue_(hdfFile,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_isSTA1")
	ITSTA5_filterList = getValue_(hdfFile,key="IceTopSTA5_12_filter")
	ITSTA3_filterList = getValue_(hdfFile,key="SDST_IceTopSTA3_12_filter")
	hlcHitList = getValue_(hdfFile,key="OfflineIceTopHLCTankPulsesTotalHit")
	evtObjList = []
	for nEvt,evtID in enumerate(eventList):
		# if STA1TriggerList[nEvt] > 0:
			# print("primary,runID,eventID,STA1,STA3,STA5",primaryList[nEvt],runList[nEvt],evtID,ITSMTTriggeredList[nEvt],STA1TriggerList[nEvt],ITSTA3_filterList[nEvt],ITSTA5_filterList[nEvt])
		thisEvt = CREvent(runList[nEvt],evtID,primaryList[nEvt],energyList[nEvt],zenithList[nEvt],coreXList[nEvt],coreYList[nEvt])
		thisEvt.addWeight(H4aWeightList[nEvt])
		# thisEvt.addTrigger(ITSMTTriggeredList[nEvt],STA1TriggerList[nEvt])
		thisEvt.addFilter(ITSTA5_filterList[nEvt],ITSTA3_filterList[nEvt])
		thisEvt.hlcHit = hlcHitList[nEvt]
		evtObjList.append(thisEvt)
	return evtObjList


def extractEventsOfficial(hdf5List):
	evtList = []
	for ihdf in hdf5List:
		evtList += extractEventsOfficial_(ihdf)
	return evtList

def containedEvents(evtList,radius):
	"""only selects the events with core with in given radius of IceTop center"""
	return [ievt for ievt in evtList if ((ievt.coreX**2 + ievt.coreY**2) <= radius**2)]


def Rdisk(energy):
	"""calculates radius of disc in m overwhich shower core spread"""
	# energy = int(np.log10(energy)*10)/10.0
	# ebin = np.log10(energy)
	ebin = energy
	# print("unit check",ebin,energy)
	#in json file "800+(int(ebin)-5)*300 + (2*(int(ebin)-5)//3)*300 + ((int(ebin)-5)//3)*300"
	#in sim file discR=$(python -c "print(800+(int($ENERGY)-5)*300 + (2*(int($ENERGY)-5)//3)*300 + ((int($ENERGY)-5)//3)*300)")
	return 800+(int(ebin)-5)*300 + (2*(int(ebin)-5)//3)*300 + ((int(ebin)-5)//3)*300

def showerArea(energy):
	'''Calculates area [km2] overwhich core of corsika shower was shifted for resimulation
	'''
	R = Rdisk(energy)
	return np.pi*R**2 * 10**(-6)

def getCircle(radius):
	theta = np.linspace(0,2*np.pi,100)
	return radius*np.cos(theta),radius*np.sin(theta)

def getVectorItem(hdfFileList,key):
	hitDuration = np.array([])
	for ihdf in hdfFileList:
		dt_hdf = pd.read_hdf(ihdf,key=key)["item"].values
		hitDuration = np.concatenate((hitDuration,dt_hdf))
	return hitDuration

def getZenith_(hdfFile):
	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
	return dataT["zenith"].values

def getCore_(hdfFile):
	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
	return dataT["x"].values,dataT["y"].values

def getCore(hdfFileList):
	xList = np.array([])
	yList = np.array([])
	for ihdf in hdfFileList:
		x,y = getCore_(ihdf)
		xList = np.concatenate((xList,x))
		yList = np.concatenate((yList,y))
	return xList,yList

def filterArray(arr,filt):
	return arr[abs(np.asarray(filt-1)) < 0.001]

def getWeight_(hdfFile):
	dataT = pd.read_hdf(hdfFile,key="H4aWeight")
	return dataT["value"].values
def getType_(hdfFile):
	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
	return dataT["type"].values
def getEnergy_(hdfFile):
	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
	return dataT["energy"].values
def getZenithTypeEnergy_(hdfFile):
	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
	return dataT["zenith"].values, dataT["type"].values, dataT["energy"].values

def getZenithTypeEnergy(hdfFileList):
	ptypeList = np.array([])
	energyList = np.array([])
	zenList = np.array([])
	for ihdf in hdfFileList:
		zen,ptype,energy = getZenithTypeEnergy_(ihdf)
		zenList = np.concatenate((zenList,zen))
		ptypeList = np.concatenate((ptypeList,ptype))
		energyList = np.concatenate((energyList,energy))
	return zenList,ptypeList,energyList

def weightCalc(hdfFileList):
	'''
	returns weights for the events in list.
	'''
	zen,ptype,energy = getZenithTypeEnergy(hdfFileList)
	weights = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zen,energy,ptype)
	adjustedWeights = []
	for iweight,ienergy in zip(weights,energy):
		if ienergy >= 10**7:
			adjustedWeights.append(iweight)
			# adjustedWeights.append(iweight*10)
		else:
			adjustedWeights.append(iweight)
	return adjustedWeights

def selectEnergyEvents(evtList,energy,tolerance):
	"tolerance in log10E, energy in eV"
	return [ievt for ievt in evtList if (np.log10(energy)-tolerance) <= np.log10(ievt.energy)+9 < (np.log10(energy)+tolerance)]


def selectTriggered(eventList,triggerType):
	"""
	select the triggered events
	"""
	if triggerType == "sta1":
		triggeredEvts =[ievt for ievt in eventList if abs(ievt.STA1Trigger-1)<0.01]
	elif triggerType == "sta3":
		triggeredEvts =[ievt for ievt in eventList if abs(ievt.ITSMTTriggered-1)<0.01]
	elif triggerType == "slc3":
		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc3Trig-1)<0.01]
	elif triggerType == "slc4":
		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc4Trig-1)<0.01]
	elif triggerType == "slc5":
		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc5Trig-1)<0.01]
	else:
		triggeredEvts = eventList
	return triggeredEvts


def triggerEfficiency(n_trig,n_total):
	# print("n_trig,n_total",n_trig,n_total)
	if n_total != 0:
		return (n_trig/n_total)
	else:
		return 0

def effectiveArea(n_trig,n_total,area):
	return triggerEfficiency(n_trig,n_total)*area


# def addDirectWeights(eventList):
# 	# nfilesP = len([ievt.runID for ievt in eventList if str(ievt.CRType) == "PPlus"])
# 	# nfilesHe = len([ievt.runID for ievt in eventList if str(ievt.CRType) == "He4Nucleus"])
# 	# nfilesO = len([ievt.runID for ievt in eventList if str(ievt.CRType) == "O16Nucleus"])
# 	# nfilesFe = len([ievt.runID for ievt in eventList if str(ievt.CRType) == "Fe56Nucleus"])
# 	nfilesP = len(set([ievt.runID for ievt in eventList if str(ievt.CRType) == str(2212)]))
# 	nfilesHe = len(set([ievt.runID for ievt in eventList if str(ievt.CRType) == str(1000020040)]))
# 	nfilesO = len(set([ievt.runID for ievt in eventList if str(ievt.CRType) == str(1000080160)]))
# 	nfilesFe = len(set([ievt.runID for ievt in eventList if str(ievt.CRType) == str(1000260560)]))
# 	print("number of files",nfilesP,nfilesHe,nfilesO,nfilesFe)
# 	zen = [ievt.zenith for ievt in eventList]
# 	energy = [ievt.energy for ievt in eventList]
# 	ptype = [ievt.CRType for ievt in eventList]
# 	weights = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zen,energy,ptype)
# 	for ievt,iweight in zip(eventList,weights):
# 		ievt.directWeight = iweight
# evtList = evtList
# # addDirectWeights(evtList)
# # zenith,ptype,energy = getZenithTypeEnergy(hdf5NullList)




# # weights = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zenith,energy,ptype)
# # weithtsOne = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zenith[-1],energy[-1],ptype[-1])
# # print("weightCompare",weights[-1],weithtsOne,len(weights))
# # adjustedWeights = []

# # for iweight,ienergy in zip(weights,energy):
# # 	if ienergy >= 10**7:
# # 		adjustedWeights.append(iweight*10)
# # 	elif ienergy < 10**6:
# # 		adjustedWeights.append(iweight*10)
# # 	else:
# # 		adjustedWeights.append(iweight)


# energyBins = 10**np.linspace(5, 8.0, 31)
# # energyBins = 10**np.linspace(5, 8.0, 7)
# energyBinslgE = np.linspace(5.0,8.9,4000)
energyBinCenter = [5.1,6.1,7.1,8.1]
# print("energy bins",energyBins)

def plotRadiusEnergy(energyBins):
	radiusList = [Rdisk(ienergy) for ienergy in energyBins]
	radiusBinValue = [Rdisk(ienergy) for ienergy in energyBinCenter]
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	ax.plot(energyBins,radiusList,"-",lw=2.5,c=qualitative_colors(3)[1],label="radius",alpha=1)
	# ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
	for x in energyBinCenter:
		ax.text(x,Rdisk(x)-200,s="{0:.0f} m".format(Rdisk(x)))
	ax.set_xlabel(r"log(Energy[GeV])", fontsize=20)
	ax.set_ylabel(r"disc radius [m]", fontsize=20)
	ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
	ax.tick_params(which='both', width=1.5)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	# ax.set_yscale('log')
	ax.grid(True,alpha=0.5)
	ax.set_ylim(0,None)
	# ax.legend(fontsize=14)
	# ax.legend(fontsize=14,ncol=2)
	ax.yaxis.set_minor_locator(MultipleLocator(100))
	ax.xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.savefig(plotFolder+"simRadius.pdf",transparent=False,bbox_inches='tight')
	plt.close()


# def plotCoreScatter_(x,y,suffix,title):
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	# ax.scatter(x,y,s=10,alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
# 	ax.set_xlabel(r"x [m]", fontsize=24)
# 	ax.set_ylabel(r"y [m]", fontsize=24)
# 	xCirc,yCirc = getCircle(800)
# 	ax.plot(xCirc,yCirc,'-',c="purple",lw=3.0,label="r = 800 m")
# 	xCirc,yCirc = getCircle(1100)
# 	ax.plot(xCirc,yCirc,'-',c="blue",lw=3.0,label="r = 1100 m")
# 	xCirc,yCirc = getCircle(1700)
# 	ax.plot(xCirc,yCirc,'-',c="orange",lw=3.0,label="r = 1700 m")
# 	# xCirc,yCirc = getCircle(2600)
# 	# ax.plot(xCirc,yCirc,'-',c="yellow",lw=3.0,label="r = 2600 m")
# 	ax.scatter(x,y,s=10,alpha=1)
# 	# ax.set_xlim(0,100)
# 	# ax.set_ylim(0,100)
# 	ax.grid(True,alpha=0.2)
# 	ax.set_title(title,fontsize=16)
# 	ax.set_aspect("equal")
# 	plt.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/coreScatter"+str(suffix)+".png",transparent=False,bbox_inches='tight')
# 	plt.close()

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

# def plotCoreScatterEnergy(evtList,energyLow,energyHigh,filtKey):
# 	x = [ievt.coreX for ievt in evtList]
# 	y = [ievt.coreY for ievt in evtList]
# 	zenithBins = [np.arcsin(np.sqrt(i)) for i in np.linspace(0.0,1.0,11)][:-3]
# 	zenithBins.append(np.deg2rad(65))
# 	print("zenith bins",[np.sin(i)**2 for i in zenithBins])
# 	for n,nZenith in enumerate(zenithBins[:-1]):
# 		xInBin = []
# 		yInBin = []
# 		for ix,iy,izen,ienergy in zip(x,y,zenList,energyList):
# 			if ienergy >= 10**energyLow and ienergy < 10**energyHigh and izen >= zenithBins[n] and izen < zenithBins[n+1]:
# 				xInBin.append(ix) 
# 				yInBin.append(iy)
# 		plotCoreScatter_(xInBin,yInBin,r"energy_{:.1f}Zen{:.1f}_filt{}".format(energyLow,np.sin(nZenith)**2,str(filtKey)),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(np.rad2deg(zenithBins[n]),np.rad2deg(zenithBins[n+1])))

# # plotCoreScatterEnergy(hdf5NullList,6.0,6.1,filtKey="IceTopSTA5_13_filter")
# # plotCoreScatterEnergy(hdf5NullList,7.0,7.1,filtKey="IceTopSTA5_13_filter")
# # plotCoreScatterEnergy(hdf5NullList,6.9,7.0,filtKey="IceTopSTA5_13_filter")
# # plotCoreScatterEnergy(hdf5NullList,6.9,7.0,filtKey="SDST_IceTopSTA3_13_filter")
# # plotCoreScatterEnergy(hdf5NullListOfficial,6.9,7.0,filtKey="IceTopSTA5_12_filter")
# # plotCoreScatterEnergy(hdf5NullList,6.0,7.0,filtKey=None)


# def plotCoreScatter_(x,y,suffix,title):
# 	fig = plt.figure(figsize=(40,25))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	# ax.scatter(x,y,s=10,alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
# 	ax.set_xlabel(r"x [m]", fontsize=24)
# 	ax.set_ylabel(r"y [m]", fontsize=24)
# 	xCirc,yCirc = getCircle(800)
# 	ax.plot(xCirc,yCirc,'-',c="purple",lw=3.0,label="r = 800 m")
# 	xCirc,yCirc = getCircle(1100)
# 	ax.plot(xCirc,yCirc,'-',c="blue",lw=3.0,label="r = 1100 m")
# 	xCirc,yCirc = getCircle(1700)
# 	ax.plot(xCirc,yCirc,'-',c="orange",lw=3.0,label="r = 1700 m")
# 	# xCirc,yCirc = getCircle(2600)
# 	# ax.plot(xCirc,yCirc,'-',c="yellow",lw=3.0,label="r = 2600 m")
# 	ax.scatter(x,y,s=10,alpha=1)
# 	# ax.set_xlim(0,100)
# 	# ax.set_ylim(0,100)
# 	ax.grid(True,alpha=0.2)
# 	ax.set_title(title,fontsize=16)
# 	ax.set_aspect("equal")
# 	plt.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/coreScatter"+str(suffix)+".png",transparent=False,bbox_inches='tight')
# 	plt.close()


# def plotScatterCore(evtList,triggerType):
# 	'''
# 	plots core distance of triggered shower
# 	'''
# 	energyBins = 10**(np.linspace(5.0,8.0,7))
# 	# evtList = selectTriggered(evtList,triggerType)
# 	for ebin, ebinStart in enumerate(energyBins[:-1]):
# 		fig = plt.figure(figsize=(8,5))
# 		gs = gridspec.GridSpec(nrows=1,ncols=1)
# 		ax = fig.add_subplot(gs[0])
# 		colorIter = iter(colorsList)
# 		lowEdge_E = energyBins[ebin]
# 		highEdge_E = energyBins[ebin+1]
# 		evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
# 		ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((energyBins[ebin])*10**9),np.log10((energyBins[ebin+1])*10**9)))
# 		for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 			# colorIter = iter(colorsCustom)
# 			lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 			highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 			evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
# 			distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
# 			x = [ievt.coreX for ievt in evtZenBin]
# 			y = [ievt.coreY for ievt in evtZenBin]
# 			if len(distanceList)>2:
# 				xbins = np.linspace(min(distanceList),max(distanceList),200)
# 			else:
# 				xbins = np.linspace(0,1000,200)
# 			ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
# 				r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),lw=2.5,alpha=1)
# 			plotCoreScatter_(x,y,r"energy_{:.1f}Zen{:.1f}".format(lowEdge_E,np.arcsin(np.sqrt(sin2ZenBins[nbin]))),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(np.arcsin(np.sqrt(sin2ZenBins[nbin])),np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))))
# 		ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 		ax.set_xlabel(r"core distance [m]", fontsize=22)
# 		ax.set_ylabel(r"count", fontsize=22)
# 		# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 		# ax.set_xscale('log')
# 		# ax.set_ylim(0,1.01)
# 		# ax.set_xlim(14,17)
# 		# ax.yaxis.set_minor_locator(MultipleLocator(100))
# 		# ax.xaxis.set_minor_locator(MultipleLocator(0.1))
# 		ax.grid(True,alpha=0.6)
# 		ax.legend(fontsize=12)
# 		plt.savefig(plotFolder+"/distanceScattTrig"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(energyBins[ebin])+9),transparent=False,bbox_inches='tight')
# 		plt.close()
# plotScatterCore(evtList,"sta3")

# def plot2dHistCoreGiven(evtList,triggerType,energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[400,500]):
# 	'''
# 	plots core distance of triggered shower
# 	'''
# 	energyBins = 10**(np.linspace(5.0,8.0,7))
# 	# evtList = selectTriggered(evtList,triggerType)
# 	fig = plt.figure(figsize=(8,8))
# 	# fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	colorIter = iter(colorsList)
# 	lowEdge_E = energyLim[0]
# 	highEdge_E = energyLim[1]
# 	evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
# 	ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((lowEdge_E)*10**9),np.log10((highEdge_E)*10**9)))
# 	# colorIter = iter(colorsCustom)
# 	lowEdge = zenLim[0]
# 	highEdge = zenLim[1]
# 	evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
# 	distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
# 	# idBug = [(ievt.runID,ievt.eventID,ievt.coreX,ievt.coreY) for ievt in evtZenBin if abs(ievt.coreX-444)<0.4 and abs(ievt.coreY-0)<1]
# 	# print("xyBug",idBug)
# 	x = [ievt.coreX for ievt in evtZenBin]
# 	y = [ievt.coreY for ievt in evtZenBin]
# 	xbins = np.linspace(min(distanceList),max(distanceList),200)
# 	# ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
# 	# 	r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(zenLim[0]*180.0/np.pi,zenLim[1]*180.0/np.pi),lw=2.5,alpha=0.4)
# 	# plotCoreScatter_(x,y,r"energy_{:.1f}Zen{:.1f}".format(lowEdge_E,zenLim[0]*180.0/np.pi),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(zenLim[0]*180.0/np.pi,zenLim[0]*180.0/np.pi))
# 	ax.hist2d(x,y,bins=100,norm=mpl.colors.LogNorm())
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 	ax.set_xlabel(r"core distance x[m]", fontsize=22)
# 	ax.set_ylabel(r"core distance y[m]", fontsize=22)
# 	ax.set_aspect("equal")
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	# ax.set_xscale('log')
# 	# ax.set_ylim(0,1.01)
# 	# ax.set_xlim(14,17)
# 	# ax.yaxis.set_minor_locator(MultipleLocator(100))
# 	# ax.xaxis.set_minor_locator(MultipleLocator(0.1))
# 	ax.grid(True,alpha=0.6)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/distance2dHist"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(lowEdge_E)+9),transparent=False,bbox_inches='tight')
# 	plt.close()
# plot2dHistCoreGiven(evtList,"sta3",energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[400,500])

# def plotScatterCoreGiven(evtList,triggerType,energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[400,500]):
# 	'''
# 	plots core distance of triggered shower
# 	'''
# 	energyBins = 10**(np.linspace(5.0,8.0,7))
# 	# evtList = selectTriggered(evtList,triggerType)
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	colorIter = iter(colorsList)
# 	lowEdge_E = energyLim[0]
# 	highEdge_E = energyLim[1]
# 	evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
# 	ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((lowEdge_E)*10**9),np.log10((highEdge_E)*10**9)))
# 	# colorIter = iter(colorsCustom)
# 	lowEdge = zenLim[0]
# 	highEdge = zenLim[1]
# 	evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
# 	distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
# 	x = [ievt.coreX for ievt in evtZenBin]
# 	y = [ievt.coreY for ievt in evtZenBin]
# 	xbins = np.linspace(min(distanceList),max(distanceList),200)
# 	ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
# 		r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(zenLim[0]*180.0/np.pi,zenLim[1]*180.0/np.pi),lw=2.5,alpha=0.4)
# 	plotCoreScatter_(x,y,r"energy_{:.1f}Zen{:.1f}".format(lowEdge_E,zenLim[0]*180.0/np.pi),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(zenLim[0]*180.0/np.pi,zenLim[0]*180.0/np.pi))
# 	# ax.hist2d(x,y,bins=100)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 	ax.set_xlabel(r"core distance [m]", fontsize=22)
# 	ax.set_ylabel(r"count", fontsize=22)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	# ax.set_xscale('log')
# 	# ax.set_ylim(0,1.01)
# 	# ax.set_xlim(14,17)
# 	# ax.yaxis.set_minor_locator(MultipleLocator(100))
# 	# ax.xaxis.set_minor_locator(MultipleLocator(0.1))
# 	ax.grid(True,alpha=0.6)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/distanceScattTrig"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(lowEdge_E)+9),transparent=False,bbox_inches='tight')
# 	plt.close()
# plotScatterCoreGiven(evtList,"sta3",energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[400,500])



# def selectTriggered(eventList,triggerType):
# 	"""
# 	select the triggered events
# 	"""
# 	if triggerType == "sta1":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.STA1Trigger-1)<0.01]
# 	elif triggerType == "sta3":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.ITSMTTriggered-1)<0.01]
# 	elif triggerType == "slc3":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc3Trig-1)<0.01]
# 	elif triggerType == "slc4":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc4Trig-1)<0.01]
# 	elif triggerType == "slc5":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc5Trig-1)<0.01]
# 	return triggeredEvts

# def addHistogram(x,ax,suffix):
#     xbins = np.linspace(min(x),max(x),200)
#     ax.hist(x,bins=xbins,histtype="step",lw=2.5,alpha=1)
#     return ax

# def plotDistanceHist(evtList,triggerType):
# 	'''
# 	plots core distance of triggered shower
# 	'''
# 	energyBins = 10**(np.linspace(5.0,8.0,7))
# 	evtList = selectTriggered(evtList,triggerType)
# 	for ebin, ebinStart in enumerate(energyBins[:-1]):
# 		fig = plt.figure(figsize=(8,5))
# 		gs = gridspec.GridSpec(nrows=1,ncols=1)
# 		ax = fig.add_subplot(gs[0])
# 		colorIter = iter(colorsList)
# 		lowEdge_E = energyBins[ebin]
# 		highEdge_E = energyBins[ebin+1]
# 		evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
# 		ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((energyBins[ebin])*10**9),np.log10((energyBins[ebin+1])*10**9)))
# 		for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 			# colorIter = iter(colorsCustom)
# 			lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 			highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 			evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
# 			distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
# 			if len(distanceList)>2:
# 				xbins = np.linspace(min(distanceList),max(distanceList),200)
# 			else:
# 				xbins = np.linspace(0,1000,200)
# 			ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
# 				r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),lw=2.5,alpha=1)
# 		ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 		ax.set_xlabel(r"core distance [m]", fontsize=22)
# 		ax.set_ylabel(r"count", fontsize=22)
# 		# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 		# ax.set_xscale('log')
# 		# ax.set_ylim(0,1.01)
# 		# ax.set_xlim(14,17)
# 		# ax.yaxis.set_minor_locator(MultipleLocator(100))
# 		# ax.xaxis.set_minor_locator(MultipleLocator(0.1))
# 		ax.grid(True,alpha=0.6)
# 		ax.legend(fontsize=12)
# 		plt.savefig(plotFolder+"/distanceHistTrig"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(energyBins[ebin])+9),transparent=False,bbox_inches='tight')
# 		plt.close()
# plotDistanceHist(evtList,triggerType="sta1")
# plotDistanceHist(evtList,triggerType="slc3")
# plotDistanceHist(evtList,triggerType="sta3")

# def plotDistanceHistGiven(evtList,triggerType,energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[400,500]):
# 	'''
# 	plots core distance of triggered shower
# 	'''
# 	energyBins = 10**(np.linspace(5.0,8.0,7))
# 	# evtList = selectTriggered(evtList,triggerType)
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	colorIter = iter(colorsList)
# 	lowEdge_E = energyLim[0]
# 	highEdge_E = energyLim[1]
# 	evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
# 	ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((lowEdge_E)*10**9),np.log10((highEdge_E)*10**9)))
# 	lowEdge = zenLim[0]
# 	highEdge = zenLim[1]
# 	evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
# 	evtZenDistBin = [ievt for ievt in evtZenBin if distanceLim[0] <= np.sqrt(ievt.coreX**2+ievt.coreY**2) < distanceLim[1]]
# 	distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenDistBin]
# 	xbins = np.linspace(min(distanceList),max(distanceList),30)
# 	evtZenDistBinTriggerd = selectTriggered(evtZenDistBin,triggerType)
# 	distanceListTrig = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenDistBinTriggerd]
# 	ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
# 		r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(zenLim[0]*180.0/np.pi,zenLim[1]*180.0/np.pi),lw=2.5,alpha=0.4)
# 	ax.hist(distanceListTrig,bins=xbins,histtype="step",color=next(colorIter),label=
# 		r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(zenLim[0]*180.0/np.pi,zenLim[1]*180.0/np.pi)+" "+str(triggerType),lw=2.5,alpha=0.4)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 	ax.set_xlabel(r"core distance [m]", fontsize=22)
# 	ax.set_ylabel(r"count", fontsize=22)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	# ax.set_xscale('log')
# 	# ax.set_ylim(0,1.01)
# 	# ax.set_xlim(14,17)
# 	# ax.yaxis.set_minor_locator(MultipleLocator(100))
# 	# ax.xaxis.set_minor_locator(MultipleLocator(0.1))
# 	ax.grid(True,alpha=0.6)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/distanceHistTrig"+str(triggerType)+"GivenEnergy{0:.1f}.pdf".format(lowEdge_E+9),transparent=False,bbox_inches='tight')
# 	plt.close()

# plotDistanceHistGiven(evtList,triggerType="sta3",energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[442,446])


# # plotCoreScatter(hdf5NullList)

# def plotZenithHist(hdfFileList):
# 	zenList,ptypeList,energyList = getZenithTypeEnergy(hdfFileList)
# 	# zenList = [np.rad2deg(i) for i in zenList]
# 	zenList = [np.sin(i)**2 for i in zenList]
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	ax.hist(zenList,label="radius",alpha=1)
# 	# ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
# 	ax.set_xlabel(r"log(Energy[GeV])", fontsize=20)
# 	ax.set_ylabel(r"disc radius [m]", fontsize=20)
# 	ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
# 	# ax.set_yscale('log')
# 	ax.grid(True,alpha=0.2)
# 	ax.set_ylim(0,None)
# 	# ax.legend(fontsize=14)
# 	# ax.legend(fontsize=14,ncol=2)
# 	plt.savefig(plotFolder+"zenHist.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # plotZenithHist(hdf5NullList)


# def plotEnergyFlux2(eventList,yscale,suffix,energyScale,triggerType):
# 	"""
# 	plots energy flux
# 	"""
# 	# triggeredEvts = [ievt.ITSTA3_filter for ievt in eventList]
# 	# print("triggers",triggeredEvts)
# 	if triggerType == "sta3":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.ITSMTTriggered-1)<0.01]
# 	elif triggerType == "sta1":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.STA1Trigger-1)<0.01]
# 	elif triggerType == "slc3":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc3Trig-1)<0.01]
# 	elif triggerType == "slc4":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc4Trig-1)<0.01]
# 	elif triggerType == "slc5":
# 		triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc5Trig-1)<0.01]
# 	energy = [ievt.energy for ievt in triggeredEvts]
# 	weights = [ievt.H4aWeight for ievt in triggeredEvts]
# 	weights_direct = [ievt.directWeight for ievt in triggeredEvts]
# 	# print("trigger",len(triggeredEvts))
# 	print("energy",energy[:5])
# 	print("weights",weights[:5])
# 	print("direct",weights_direct[:5])
# 	energy = np.log10(energy)+9.0
# 	# print("lengths",len(weights),weights,len(energy))
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	hitBins = np.linspace(14.0,17.0,31)
# 	hist,binEdge = np.histogram(energy,hitBins,weights=[w for w in weights])
# 	# binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
# 	print("sum of simulated rate",sum(hist))
# 	ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
# 	binCenter = (binEdge[:-1]+binEdge[1:])/2.0
# 	if str(energyScale) == "0.0":
# 		H = [h for h in hist]
# 		ax.set_ylim(10**-4,20)
# 	else:
# 		H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
# 	# ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
# 	ax.step(binCenter,H,"-",where="mid",lw=2.5,label="after weighting",alpha=1)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
# 	ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
# 	ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20)
# 	ax.set_yscale(yscale)	
# 	# ax.set_ylim(None,10**5)
# 	# ax.set_xscale('log')
# 	ax.grid(True,alpha=0.2)
# 	# ax.set_title(key,fontsize=16)
# 	ax.legend(fontsize=10,ncol=3,loc="lower center")
# 	plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+".pdf",transparent=False,bbox_inches='tight')
# 	plt.close()
# # plotEnergyFlux(energy,weights,adjustedWeights,"log","flux")
# # plotEnergyFlux(evtList,"linear","fluxlinear",1.8,triggerType="sta3")
# # plotEnergyFlux2(evtList,"linear","fluxlinear",0.0,triggerType="sta3")
# # plotEnergyFlux2(evtList,"linear","fluxlinear",0.0,triggerType="sta1")
# # plotEnergyFlux2(evtList,"log","fluxLog",0.0,triggerType="sta1")
# # plotEnergyFlux2(evtList,"log","fluxLog",0.0,triggerType="sta3")
# # plotEnergyFlux2(evtList,"log","fluxLog",1.8,triggerType="sta1")
# # plotEnergyFlux2(evtList,"log","fluxLog",1.8,triggerType="sta3")


# def plotSteps(triggeredEvts,ax,legendLabel,energyScale):
# 	hitBins = np.linspace(14.0,17.0,31)
# 	energy = [ievt.energy for ievt in triggeredEvts]
# 	weights = [ievt.H4aWeight for ievt in triggeredEvts]
# 	weights_direct = [ievt.directWeight for ievt in triggeredEvts]
# 	energy = np.log10(energy)+9.0
# 	hist,binEdge = np.histogram(energy,hitBins,weights=[w for w in weights])
# 	# binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
# 	print("sum of simulated rate",sum(hist))
# 	# if "total" in legendLabel:
# 		# ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
# 	# ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
# 		# ax.set_ylim(10**-4,20)
# 	binCenter = (binEdge[:-1]+binEdge[1:])/2.0
# 	if str(energyScale) == "0.0":
# 		H = [h for h in hist]		
# 	else:
# 		H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
# 	# ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
# 	ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.1f} Hz".format(sum(hist)),alpha=1)
# 	return ax
def plotRatioSteps(triggeredEvts1,triggeredEvts2,ax,legendLabel,ncolor):
	hitBins = np.linspace(14.0,17.0,31)
	energy1 = [ievt.energy for ievt in triggeredEvts1]
	weights1 = [ievt.H4aWeight for ievt in triggeredEvts1]
	# weights_direct1 = [ievt.directWeight for ievt in triggeredEvts1]
	energy1 = np.log10(energy1)+9.0
	hist1,binEdge1 = np.histogram(energy1,hitBins,weights=[w for w in weights1])
	energy2 = [ievt.energy for ievt in triggeredEvts2]
	weights2 = [ievt.H4aWeight for ievt in triggeredEvts2]
	# weights_direct2 = [ievt.directWeight for ievt in triggeredEvts2]
	energy2 = np.log10(energy2)+9.0
	hist2,binEdge2 = np.histogram(energy2,hitBins,weights=[w for w in weights2])
	print("edge compares",binEdge1,binEdge2)
	histRatio = [0 if abs(j-0) < 0.000000001 else i/j for i,j in zip(hist2,hist1)]
	# binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
	binCenter = (binEdge1[:-1]+binEdge1[1:])/2.0
	H = [h for h in histRatio]		
	ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel,color=ncolor,alpha=1)
	return ax


# def plotEnergyFlux(eventList,yscale,suffix,energyScale,triggerType,plotTotal=True,sepZenBins=True):
# 	"""
# 	plots energy flux
# 	"""
# 	# triggeredEvts = [ievt.ITSTA3_filter for ievt in eventList]
# 	# print("triggers",triggeredEvts)
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	hitBins = np.linspace(14.0,17.0,31)
# 	if sepZenBins == False:
# 		if triggerType == "sta3":
# 			triggeredEvts =[ievt for ievt in eventList if abs(ievt.ITSMTTriggered-1)<0.01]
# 		elif triggerType == "sta1":
# 			triggeredEvts =[ievt for ievt in eventList if abs(ievt.STA1Trigger-1)<0.01]
# 		elif triggerType == "slc3":
# 			triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc3Trig-1)<0.01]
# 		elif triggerType == "slc4":
# 			triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc4Trig-1)<0.01]
# 		elif triggerType == "slc5":
# 			triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc5Trig-1)<0.01]
# 		print("debug no of triggered events",len(triggeredEvts),triggerType)
# 		if plotTotal==True:
# 			ax = plotSteps(eventList,ax,"{} events".format("total"),energyScale)
# 		ax = plotSteps(triggeredEvts,ax,"{} trigger".format(triggerType),energyScale)
# 	elif sepZenBins == True:
# 		for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 			lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 			highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 			evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
# 			energyList = []
# 			neventList = []
# 			if triggerType == "sta3":
# 				triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.ITSMTTriggered-1)<0.01]
# 			elif triggerType == "sta1":
# 				triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.STA1Trigger-1)<0.01]
# 			elif triggerType == "slc3":
# 				triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.slc3Trig-1)<0.01]
# 			elif triggerType == "slc4":
# 				triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.slc4Trig-1)<0.01]
# 			elif triggerType == "slc5":
# 				triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.slc5Trig-1)<0.01]
# 			else:
# 				print("unkbown trigger type",triggerType)
# 			if plotTotal==True:
# 				ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
# 			ax = plotSteps(triggeredEvts,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),energyScale)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
# 	ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
# 	ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20)
# 	ax.set_yscale(yscale)	
# 	# ax.set_ylim(None,10**5)
# 	# ax.set_xscale('log')
# 	ax.grid(True,alpha=0.2)
# 	# ax.set_title(key,fontsize=16)
# 	ax.legend(fontsize=10,ncol=3,loc="lower center")
# 	plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"Total"+str(plotTotal)+"ZenBin"+str(sepZenBins)+".pdf",transparent=False,bbox_inches='tight')
# 	plt.close()
# # plotEnergyFlux(energy,weights,adjustedWeights,"log","flux")
# # plotEnergyFlux(evtList,"linear","fluxlinear",1.8,triggerType="sta3")
# # plotEnergyFlux(evtList,"linear","fluxlinear",0.0,triggerType="sta3")
# # plotEnergyFlux(evtList,"linear","fluxlinear",0.0,triggerType="sta1")
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="sta1",plotTotal=False,sepZenBins=True)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="sta1",plotTotal=True,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="sta1",plotTotal=False,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="sta3",plotTotal=False,sepZenBins=True)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="sta3",plotTotal=True,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="sta3",plotTotal=False,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",1.8,triggerType="sta3",plotTotal=False,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",1.8,triggerType="sta1",plotTotal=False,sepZenBins=True)
# # plotEnergyFlux(evtList,"log","fluxLog",1.8,triggerType="sta1",plotTotal=False,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",1.8,triggerType="sta3",plotTotal=True,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc3",plotTotal=False,sepZenBins=True)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc3",plotTotal=True,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc3",plotTotal=False,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",1.8,triggerType="slc3",plotTotal=False,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc4",plotTotal=False,sepZenBins=True)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc4",plotTotal=True,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc4",plotTotal=False,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc5",plotTotal=False,sepZenBins=True)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc5",plotTotal=True,sepZenBins=False)
# # plotEnergyFlux(evtList,"log","fluxLog",0.0,triggerType="slc5",plotTotal=False,sepZenBins=False)

# HLCStationDuration = getVectorItem(hdf5NullList,"OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t")

# evtListOfficial = extractEventsOfficial(hdf5NullListOfficial)

# def plotEnergyFluxRatio(eventList,yscale,suffix,triggerType,sepZenBins=True):
# 	"""
# 	plots energy flux
# 	"""
# 	# triggeredEvts = [ievt.ITSTA3_filter for ievt in eventList]
# 	# print("triggers",triggeredEvts)
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	hitBins = np.linspace(14.0,17.0,31)
# 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 		lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 		highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 		evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
# 		energyList = []
# 		neventList = []
# 		triggeredEvts1 =[ievt for ievt in evtZenBin if abs(ievt.ITSMTTriggered-1)<0.01]
# 		if triggerType == "sta1":
# 			triggeredEvts2 =[ievt for ievt in evtZenBin if abs(ievt.STA1Trigger-1)<0.01]
# 		elif triggerType == "slc3":
# 			triggeredEvts2 =[ievt for ievt in evtZenBin if abs(ievt.slc3Trig-1)<0.01]
# 		elif triggerType == "slc4":
# 			triggeredEvts2 =[ievt for ievt in evtZenBin if abs(ievt.slc4Trig-1)<0.01]
# 		elif triggerType == "slc5":
# 			triggeredEvts2 =[ievt for ievt in evtZenBin if abs(ievt.slc5Trig-1)<0.01]
# 		else:
# 			print("unknown trigger type",triggerType)
# 		ax = plotRatioSteps(triggeredEvts1,triggeredEvts2,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi))
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
# 	ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
# 	ax.set_ylabel(r"rate $\frac{{{0}}}{{{1}}}$".format(triggerType,"sta3"), fontsize=20)
# 	ax.set_yscale(yscale)	
# 	ax.set_ylim(0.1,None)
# 	# ax.set_xscale('log')
# 	ax.grid(True,alpha=0.2)
# 	# ax.set_title(key,fontsize=16)
# 	ax.legend(fontsize=10,ncol=3)
# 	plt.savefig(plotFolder+"/energySpecRatioTrig"+str(triggerType)+str(suffix)+"ZenBin"+str(sepZenBins)+".pdf",transparent=False,bbox_inches='tight')
# 	plt.close()
# # plotEnergyFluxRatio(evtList,"log","fluxLog",triggerType="slc3",sepZenBins=True)
# # plotEnergyFluxRatio(evtList,"log","fluxLog",triggerType="sta1",sepZenBins=True)


# def plotEnergySpectraCompare(evtList,evtListOfficial,filtKey,filtKeyOfficial,yscale,suffix,energyScale):
# 	if filtKey == "IceTopSTA5_13_filter":
# 		evtListTrig = [ievt for ievt in evtList if abs(ievt.ITSTA5_filter-1)<0.01]
# 	elif filtKey == "SDST_IceTopSTA3_13_filter":
# 		evtListTrig = [ievt for ievt in evtList if abs(ievt.ITSTA3_filter-1)<0.01]
# 	if filtKeyOfficial == "IceTopSTA5_12_filter":
# 		evtListTrigOfficial = [ievt for ievt in evtListOfficial if abs(ievt.ITSTA5_filter-1)<0.01]
# 	elif filtKeyOfficial == "SDST_IceTopSTA3_12_filter":
# 		evtListTrigOfficial = [ievt for ievt in evtListOfficial if abs(ievt.ITSTA3_filter-1)<0.01]
# 	else:
# 		print("filtkey",filtKeyOfficial)
# 	energy = np.log10([ievt.energy for ievt in evtListTrig])
# 	weights = [ievt.H4aWeight for ievt in evtListTrig]
# 	energyOfficial = np.log10([ievt.energy for ievt in evtListTrigOfficial])
# 	weightsOfficial = [ievt.H4aWeight for ievt in evtListTrigOfficial]
# 	# print("lengths",len(weights),weights,len(energy))
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	hitBins = np.linspace(5.0,8.0,31)
# 	hist,binEdge = np.histogram(energy,hitBins,weights=weights)
# 	histOfficial,binEdge = np.histogram(energyOfficial,hitBins,weights=weightsOfficial)
# 	print("sum of official rate",sum(histOfficial))
# 	print("sum of simulated rate",sum(hist))
# 	ax.text(0.05, 0.95, r"{} official rate:{:.1f} Hz, sim. rate:{:.1f} Hz".format(filtKeyOfficial,sum(histOfficial),sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
# 	binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
# 	H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
# 	HOfficial = [h*(10**E)**energyScale for h,E in zip(histOfficial,binCenter)]
# 	######################################################################################
# 	# ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
# 	# ax.hist(energyOfficial,bins=hitBins,histtype="step",lw=2.5,label=r"after weighting",alpha=1)
# 	####################################################################################
# 	ax.step(binCenter,H,"-",where="mid",lw=2.5,label="my simulation",alpha=0.6)
# 	ax.step(binCenter,HOfficial,"-",where="mid",lw=2.5,label="official simulation",alpha=0.6)
# 	####################################################################################
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
# 	ax.set_xlabel(r"log10(energy [GeV])", fontsize=20)
# 	ax.set_ylabel(r"$E^{0:.1f}$ rate [Hz]".format(energyScale), fontsize=20)
# 	ax.set_yscale(yscale)
# 	# ax.set_ylim(10**10,10**11.5)
# 	# ax.set_ylim(10**6,10**12)
# 	# ax.set_xlim(6.5,7.5)
# 	# ax.set_xscale('log')
# 	ax.grid(True,alpha=0.2)
# 	# ax.set_title(key,fontsize=16)
# 	ax.legend(fontsize=10,ncol=3,loc="lower center")
# 	plt.savefig(plotFolder+"/energySpec"+str(suffix)+"Compare.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # plotEnergySpectraCompare(evtList,evtListOfficial,"IceTopSTA5_13_filter","IceTopSTA5_12_filter","log","fluxLogSTA5",energyScale=1.8)
# # plotEnergySpectraCompare(evtList,evtListOfficial,"IceTopSTA5_13_filter","IceTopSTA5_12_filter","log","fluxLogSTA5",energyScale=1.8)
# # plotEnergySpectraCompare(evtList,evtListOfficial,"SDST_IceTopSTA3_13_filter","SDST_IceTopSTA3_12_filter","log","fluxLogSTA3",energyScale=1.8)


# def delta_t_hist(eventList,histBins,suffix):
# 	if "HLC" in suffix:
# 		deltaTList = [ievt.deltaTHLC for ievt in eventList]
# 	elif "SLC" in suffix:
# 		deltaTList = [ievt.deltaT3SLC for ievt in eventList]
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	hitBins = np.linspace(min(deltaTList),max(deltaTList),200)
# 	ax.hist(deltaTList,bins=histBins,histtype="step",lw=2.5,alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
# 	ax.set_xlabel(r"$\Delta_t$ [ns]", fontsize=20)
# 	ax.set_ylabel(r"count", fontsize=20)
# 	ax.set_yscale('log')
# 	ax.set_xlim(0,10000)
# 	ax.grid(True,alpha=0.2)
# 	ax.set_title(suffix,fontsize=12)
# 	# ax.legend(fontsize=20)
# 	plt.savefig(plotFolder+"/"+str(suffix)+"delta_t.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()
# # delta_t_hist(evtList,"HLCVEM")
# # delta_t_hist(evtList,np.linspace(0,10000,2000),"SLCVEM")

# def delta_t_hist_zen_bins(eventList,histBins,suffix):
# 	if "HLC" in suffix:
# 		deltaTList = [ievt.deltaTHLC for ievt in eventList]
# 	elif "SLC" in suffix:
# 		deltaTList = [ievt.deltaT3SLC for ievt in eventList]
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	# hitBins = np.linspace(min(deltaTList),max(deltaTList),200)
# 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 		lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 		highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 		thisEvtBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
# 		if "HLC" in suffix:
# 			deltaTList = [ievt.deltaTHLC for ievt in thisEvtBin]
# 		elif "SLC" in suffix:
# 			deltaTList = [ievt.deltaT3SLC for ievt in thisEvtBin]
# 		ax.hist(deltaTList,bins=histBins,
# 			histtype="step",label=r"$\theta$ = {0:.0f}$^{{\circ}}$-{1:.0f}$^{{\circ}}$, {2:d} evts".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,len(deltaTList)),color=next(colorsIter),lw=2.5,alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
# 	ax.set_xlabel(r"$\Delta_t$ [ns]", fontsize=20)
# 	ax.set_ylabel(r"count", fontsize=20)
# 	ax.set_yscale('log')
# 	ax.set_xlim(0,3000)
# 	ax.grid(True,alpha=0.2)
# 	ax.set_title(suffix,fontsize=12)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/"+str(suffix)+"ZenBinsdelta_t.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()
# # delta_t_hist(evtList,"HLCVEM")
# # delta_t_hist_zen_bins(evtList,np.linspace(0,3000,600),"SLCVEM")


# # delta_t_hist_zen_bins(hdf5NullList,[0,10,20,30,40,50,65.1],"HLCVEM","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t")
# # delta_t_hist_zen_bins(hdf5NullList,[0,10,20,30,40,50,65.1],"SLCVEM","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_delta_t")

# def getEventsZenith_(hdfFile,zenLim):
# 	'''get events belonging to given zenith bin:
# 	zenith angles in degree
# 	'''	
# 	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
# 	events_MC = df_MCPrimary["Event"].values
# 	zeniths_MC = df_MCPrimary["zenith"].values
# 	selEvents = []
# 	for ievent,izenith in zip(events_MC,zeniths_MC):
# 		if np.rad2deg(izenith) >= zenLim[0] and np.rad2deg(izenith) < zenLim[1]:
# 			selEvents.append(ievent)
# 	return list(set(selEvents))

# def getTriggeredEvents_(hdfFile,zenLim,energyLim,weighting):
# 	"""
# 	returns total events sta1 and sta3 trigered events
# 	"""
# 	STA1Trigger_df = pd.read_hdf(hdfFile,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_isSTA1")	
# 	STA3Trigger_df = pd.read_hdf(hdfFile,key="ITSMTTriggered")
# 	# selEvents = getEventsZenith_(hdfFile,zenLim)
# 	selEvents = getEventsZenithEnergy_(hdfFile,zenLim,energyLim)
# 	eventList = STA1Trigger_df["Event"].values
# 	totalEvts = np.ones(len(eventList))
# 	if weighting == True:
# 		# weights,_ = weightCalc(hdfFile)
# 		weights = getValue_(hdfFile,"H4aWeight")
# 		totalEvts = np.multiply(totalEvts,weights)
# 		STA1Triggers = np.multiply(STA1Trigger_df["value"].values,weights)
# 		STA3Triggers = np.multiply(STA3Trigger_df["value"].values,weights)
# 	else:
# 		STA1Triggers = STA1Trigger_df["value"].values
# 		STA3Triggers = STA3Trigger_df["value"].values
# 	STA1Triggers_select = []
# 	STA3Triggers_select = []
# 	totalEvts_select = []
# 	for ievt,ista1,ista3,isEvt in zip(eventList,STA1Triggers,STA3Triggers,totalEvts):
# 		if ievt in selEvents:
# 			STA1Triggers_select.append(ista1)
# 			STA3Triggers_select.append(ista3)
# 			totalEvts_select.append(isEvt)
# 	return np.sum(totalEvts),np.sum(STA1Triggers_select),np.sum(STA3Triggers_select)
	
# def getTriggeredEvents(hdfFileList,zenLim,energyLim,weighting):
# 	"""
# 	sums total sta1 and sta3 events from all hdfFileList files
# 	"""
# 	totalEvts = 0
# 	sta1_evts = 0
# 	sta3_evts = 0
# 	for ihdf in hdfFileList:
# 		n_total,n_sta1,n_sta3 = getTriggeredEvents_(ihdf,zenLim,energyLim,weighting)
# 		totalEvts += n_total
# 		sta1_evts += n_sta1
# 		sta3_evts += n_sta3
# 	return totalEvts,sta1_evts,sta3_evts

# def getEventsZenithEnergy_(hdfFile,zenLim,energyLim):
# 	'''get events belonging to given zenith bin:
# 	zenith angles in degree
# 	'''	
# 	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
# 	events_MC = df_MCPrimary["Event"].values
# 	zeniths_MC = df_MCPrimary["zenith"].values
# 	energy_MC = df_MCPrimary["energy"].values
# 	selEvents = []
# 	for ievent,izenith,ienergy in zip(events_MC,zeniths_MC,energy_MC):
# 		if np.rad2deg(izenith) >= zenLim[0] and np.rad2deg(izenith) < zenLim[1]:
# 			if ienergy >= energyLim[0] and ienergy < energyLim[1]:
# 				selEvents.append(ievent)
# 	return list(set(selEvents))


# def trigZen(hdfFileList,zenBins,fraction,weighting,suffix):
# 	"""
# 	plots number of triggered events in different zenith angle bins
# 	"""
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	angleBins = []
# 	sta1_list = []
# 	sta3_list = []
# 	totalEvts_list = []
# 	sta1Frac_list = []
# 	sta3Frac_list = []
# 	for n,bin_start in enumerate(zenBins[:-1]):
# 		lowEdge = zenBins[n]
# 		highEdge = zenBins[n+1]
# 		totalEvts,sta1_evts,sta3_evts = getTriggeredEvents(hdfFileList,[lowEdge,highEdge],[5.0,8.0],weighting)
# 		angleBins.append((lowEdge+highEdge)/2.0)
# 		sta1_list.append(sta1_evts)
# 		sta3_list.append(sta3_evts)
# 		totalEvts_list.append(totalEvts)
# 	if fraction == True:
# 		sta1Frac_list = np.asarray(sta1_list)/np.asarray(totalEvts_list)
# 		sta3Frac_list = np.asarray(sta3_list)/np.asarray(totalEvts_list)
# 		ax.plot(angleBins,sta1Frac_list,"o-",c=next(colorsIter),label="STA1",alpha=1)
# 		ax.plot(angleBins,sta3Frac_list,"o-",c=next(colorsIter),label="STA3",alpha=1)
# 		ax.set_ylabel(r"fraction", fontsize=20)
# 	elif fraction == False:
# 	    ax.plot(angleBins,sta1_list,"o-",c=next(colorsIter),label="STA1",alpha=1)
# 	    ax.plot(angleBins,sta3_list,"o-",c=next(colorsIter),label="STA3",alpha=1)
# 	    # ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
# 	    ax.set_ylabel(r"rate [Hz]", fontsize=20)
# 	    for iangle,ista1,ista3,itot in zip(angleBins,sta1_list,sta3_list,totalEvts_list):
# 	    	# ax.text(iangle,ista1+1,s="{0:.3f}".format(ista1/itot))
# 	    	# ax.text(iangle,ista3+1,s="{0:.3f}".format(ista3/itot))
# 	    	# ax.text(iangle,itot-200,s="{0:.0f}".format(itot))
# 	    	ax.text(iangle,ista1,s="{0:.4f}%".format(ista1/itot*100))
# 	    	ax.text(iangle,ista3,s="{0:.4f}%".format(ista3/itot*100))
# 	    	# ax.text(iangle,itot,s="{0:.0f}".format(itot))
# 	ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
# 	ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
# 	# ax.set_yscale('log')
# 	ax.set_title("total rate: STA1 = {:.2f} Hz, STA3 = {:.2f} Hz".format(np.sum(sta1_list),np.sum(sta3_list)),fontsize=12)
# 	ax.set_xlim(0,65)
# 	ax.grid(True,alpha=0.2)
# 	ax.legend(fontsize=14)
# 	# ax.legend(fontsize=14,ncol=2)
# 	plt.savefig(plotFolder+"/"+str(suffix)+"Frac"+str(fraction)+"Wts"+str(weighting)+"trigger_zen.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()


# # trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,False,"HLCVEM")
# # trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],False,False,"HLCVEM")
# # trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,True,"HLCVEM")
# # trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],False,True,"HLCVEM")





# def plotTrigEfficiency(evtList,energyBins,weighting,triggerType,containment):
# 	'''
# 	plots trigger efficiency in different zenith bins
# 	'''
# 	if containment == True:
# 		evtList = containedEvents(evtList,640)
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	# colorIter = iter(colorsCustom)
# 	colorIter = iter(colorsList)
# 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 		lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 		highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 		evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
# 		energyList = []
# 		efficiencyList = []
# 		for ebin, ebinStart in enumerate(energyBins[:-1]):
# 			lowEdge_E = energyBins[ebin]
# 			highEdge_E = energyBins[ebin+1]
# 			evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
# 			# totalEvts = len(evtEBin)
# 			weights = [ievt.H4aWeight for ievt in evtEBin]
# 			totalEvts = len(evtEBin)
# 			# sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
# 			if triggerType == "sta3":
# 				triggerList = [ievt.ITSMTTriggered for ievt in evtEBin]
# 			elif triggerType == "sta1":
# 				triggerList = [ievt.STA1Trigger for ievt in evtEBin]
# 			elif triggerType == "slc3":
# 				triggerList = [ievt.slc3Trig for ievt in evtEBin]
# 			elif triggerType == "slc4":
# 				triggerList = [ievt.slc4Trig for ievt in evtEBin]
# 			elif triggerType == "slc5":
# 				triggerList = [ievt.slc5Trig for ievt in evtEBin]
# 			trigEff = triggerEfficiency(sum(triggerList),totalEvts)
# 			efficiencyList.append(trigEff)
# 			energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))		
# 		ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 	ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
# 	ax.set_ylabel(r"trigger efficiency", fontsize=22)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	# ax.set_xscale('log')
# 	ax.set_ylim(0,1.01)
# 	ax.set_xlim(14,17)
# 	# ax.yaxis.set_minor_locator(MultipleLocator(100))
# 	ax.xaxis.set_minor_locator(MultipleLocator(0.1))
# 	ax.grid(True,alpha=0.5)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"Efficiency.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # plotTrigEfficiency(evtList,energyBins,weighting=False)
# # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="sta1",containment=False)
# # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="sta3",containment=False)
# # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="slc3",containment=False)
# # # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="slc4",containment=False)
# # # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="slc5",containment=False)
# # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="sta1",containment=True)
# # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="sta3",containment=True)
# # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="slc3",containment=True)
# # # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="slc4",containment=True)
# # # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType="slc5",containment=True)

# def plotEffectiveArea(evtList,energyBins,weighting):
# 	'''
# 	plots trigger efficiency in different zenith bins
# 	'''
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	# colorIter = iter(colorsCustom)
# 	colorIter = iter(colorsList)
# 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 		lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 		highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 		evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
# 		energyList = []
# 		sta3_efficiencyList = []
# 		sta1_efficiencyList = []
# 		for ebin, ebinStart in enumerate(energyBins[:-1]):
# 			lowEdge_E = energyBins[ebin]
# 			highEdge_E = energyBins[ebin+1]
# 			evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
# 			# totalEvts = len(evtEBin)
# 			weights = [ievt.H4aWeight for ievt in evtEBin]
# 			# totalEvts = np.sum(weights)
# 			totalEvts = len(evtEBin)
# 			# sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
# 			sta3 = [ievt.ITSMTTriggered for ievt in evtEBin]
# 			# sta3 = [ievt.ITSTA5_filter for ievt in evtEBin]
# 			sta1 = [ievt.STA1Trigger for ievt in evtEBin]
# 			sta3_trigEff = effectiveArea(sum(sta3),totalEvts,showerArea(np.log10(lowEdge_E)))
# 			sta1_trigEff = triggerEfficiency(sum(sta1),totalEvts)
# 			sta3_efficiencyList.append(sta3_trigEff)
# 			sta1_efficiencyList.append(sta1_trigEff)
# 			energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))	
# 		ax.plot(energyList,sta3_efficiencyList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 	ax.set_xlabel(r"log10 (energy [eV])", fontsize=22)
# 	ax.set_ylabel(r"effective area [$km^{2}$]", fontsize=22)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	ax.set_xscale('log')
# 	# ax.set_ylim(0,4)
# 	ax.grid(True,alpha=0.2)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/effectiveArea.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # plotEffectiveArea(evtList,energyBins,weighting=False)
# # plotEffectiveArea(evtList_contained,energyBins,weighting=False)


# def plotFluxEnergyZenith(hdfFileList,energyBins,weighting):
# 	'''
# 	plots trigger efficiency in different zenith bins
# 	'''
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	# colorIter = iter(colorsCustom)
# 	colorIter = iter(colorsList)
# 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 		lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 		highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 		evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
# 		energyList = []
# 		sta3_fluxList = []
# 		sta1_fluxList = []
# 		totalEvts_list = []
# 		for ebin, ebinStart in enumerate(energyBins[:-1]):
# 			lowEdge_E = energyBins[ebin]
# 			highEdge_E = energyBins[ebin+1]
# 			evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
# 			# totalEvts = len(evtEBin)
# 			# weights = [ievt.H4aWeight for ievt in evtEBin]
# 			weights = [ievt.directWeight for ievt in evtEBin]
# 			totalEvts = np.sum(weights)
# 			sta3_flux = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
# 			sta1_flux = [ievt.STA1Trigger*ievt.H4aWeight for ievt in evtEBin]
# 			# sta3_fluxList.append(sum(sta3_flux))
# 			sta3_fluxList.append(sum(sta3_flux)/totalEvts*100.0)
# 			totalEvts_list.append(totalEvts)
# 			energyList.append((lowEdge_E+highEdge_E)/2.0)
# 			# print("zen energy frac",binStart,np.log10(ebinStart),sum(sta3_flux)/totalEvts*100.0,sum(sta1_flux)/totalEvts*100.0)
# 			# if sum(sta3_flux) > 10**-3:
# 			# 	ax.text((lowEdge_E+highEdge_E)/2.0,sum(sta3_flux),s="{0:.2f}".format(sum(sta3_flux)/totalEvts*100.0),fontsize=8)		
# 		refLine=ax.plot(energyList,sta3_fluxList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# 		ax.plot(energyList,totalEvts_list,".",ls='-',lw = 2.5,c=refLine[0].get_color(),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 	ax.set_xlabel(r"energy [GeV]", fontsize=22)
# 	# ax.set_ylabel(r"rate [Hz]", fontsize=22))
# 	ax.set_ylabel(r"percentage", fontsize=22)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	ax.set_xscale('log')
# 	# ax.set_yscale('log')
# 	ax.grid(True,alpha=0.2)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/fluxEnergyZenith.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # plotFluxEnergyZenith(hdf5NullList,energyBins,weighting=False)


# def plotFluxEnergyZenithScaled(hdfFileList,energyBins,weighting):
# 	'''
# 	plots trigger efficiency in different zenith bins
# 	'''
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	# colorIter = iter(colorsCustom)
# 	colorIter = iter(colorsList)
# 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 		lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 		highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 		evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
# 		energyList = []
# 		totalEvts_list = []
# 		sta3_fluxList = []
# 		sta1_fluxList = []
# 		for ebin, ebinStart in enumerate(energyBins[:-1]):
# 			lowEdge_E = energyBins[ebin]
# 			highEdge_E = energyBins[ebin+1]
# 			evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
# 			# totalEvts = len(evtEBin)
# 			weights = [ievt.H4aWeight for ievt in evtEBin]
# 			totalEvts = np.sum(weights)
# 			sta3_flux = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
# 			sta3_fluxList.append(sum(sta3_flux)*((lowEdge_E+highEdge_E)/2.0)**1.8)
# 			totalEvts_list.append(totalEvts*((lowEdge_E+highEdge_E)/2.0)**1.8)
# 			energyList.append((lowEdge_E+highEdge_E)/2.0)		
# 		refLine=ax.plot(energyList,sta3_fluxList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# 		ax.plot(energyList,totalEvts_list,".",ls='-',lw = 2.5,c=refLine[0].get_color(),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 	ax.set_xlabel(r"energy [GeV]", fontsize=22)
# 	ax.set_ylabel(r"$E^{1.8}$ rate [Hz]", fontsize=22)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	ax.set_xscale('log')
# 	ax.set_yscale('log')
# 	ax.grid(True,alpha=0.2)
# 	# ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/fluxEnergyZenithScaled.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # plotFluxEnergyZenithScaled(hdf5NullList,energyBins,weighting=False)



# # def plotTrigEfficiency(hdfFileList,energyBins,weighting):
# # 	'''
# # 	plots trigger efficiency in different zenith bins
# # 	'''
# # 	fig = plt.figure(figsize=(8,5))
# # 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# # 	ax = fig.add_subplot(gs[0])
# # 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# # 		lowEdge = sin2ZenBins[nbin]
# # 		highEdge = sin2ZenBins[nbin+1]
# # 		energyList = []
# # 		sta3_efficiencyList = []
# # 		sta1_efficiencyList = []
# # 		for ebin, ebinStart in enumerate(energyBins[:-1]):
# # 			lowEdge_E = energyBins[ebin]
# # 			highEdge_E = energyBins[ebin+1]
# # 			totalEvts,sta1_evts,sta3_evts = getTriggeredEvents(hdfFileList,[lowEdge,highEdge],[lowEdge_E,highEdge_E],weighting)
# # 			sta3_trigEff = triggerEfficiency(sta3_evts,totalEvts)
# # 			sta1_trigEff = triggerEfficiency(sta1_evts,totalEvts)
# # 			sta3_efficiencyList.append(sta3_trigEff)
# # 			sta1_efficiencyList.append(sta1_trigEff)
# # 			energyList.append((lowEdge_E+highEdge_E)/2.0)		
# # 		ax.plot(energyList,sta3_efficiencyList,".",ls='-',lw = 2.5,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# # 	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
# # 	ax.set_xlabel(r"energy [GeV]", fontsize=24)
# # 	ax.set_ylabel(r"trigger efficiency", fontsize=24)
# # 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# # 	ax.set_xscale('log')
# # 	ax.grid(True,alpha=0.2)
# # 	ax.legend(fontsize=12)
# # 	plt.savefig(plotFolder+"/trigEfficiency.pdf",transparent=False,bbox_inches='tight')
# # 	plt.close()

# # plotTrigEfficiency(hdf5NullList,energyBins,weighting=False)

# def plotEffectiveArea(evtList,energyBins,weighting,triggerType,containment):
# 	'''
# 	plots trigger efficiency in different zenith bins
# 	'''
# 	if containment == True:
# 		evtList = containedEvents(evtList,640)
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	# colorIter = iter(colorsCustom)
# 	colorIter = iter(colorsList)
# 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 		lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
# 		highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
# 		evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
# 		energyList = []
# 		efficiencyList = []
# 		for ebin, ebinStart in enumerate(energyBins[:-1]):
# 			lowEdge_E = energyBins[ebin]
# 			highEdge_E = energyBins[ebin+1]
# 			evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
# 			# totalEvts = len(evtEBin)
# 			weights = [ievt.H4aWeight for ievt in evtEBin]
# 			totalEvts = len(evtEBin)
# 			# sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
# 			if triggerType == "sta3":
# 				triggerList = [ievt.ITSMTTriggered for ievt in evtEBin]
# 			elif triggerType == "sta1":
# 				triggerList = [ievt.STA1Trigger for ievt in evtEBin]
# 			elif triggerType == "slc3":
# 				triggerList = [ievt.slc3Trig for ievt in evtEBin]
# 			elif triggerType == "slc4":
# 				triggerList = [ievt.slc4Trig for ievt in evtEBin]
# 			elif triggerType == "slc5":
# 				triggerList = [ievt.slc5Trig for ievt in evtEBin]
# 			effArea = effectiveArea(sum(triggerList),totalEvts,showerArea(np.log10(lowEdge_E)))
# 			efficiencyList.append(effArea)
# 			energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))		
# 		ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
# 	ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
# 	ax.set_ylabel(r"effective area [$km^{2}$]", fontsize=22)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	# ax.set_xscale('log')
# 	ax.set_ylim(0,5)
# 	ax.set_xlim(14,17)
# 	# ax.yaxis.set_minor_locator(MultipleLocator(100))
# 	ax.xaxis.set_minor_locator(MultipleLocator(0.1))
# 	ax.grid(True,alpha=0.5)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"EffectiveArea.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="sta1",containment=False)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="sta3",containment=False)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc3",containment=False)
# # # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc4",containment=False)
# # # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc5",containment=False)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="sta1",containment=True)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="sta3",containment=True)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc3",containment=True)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc4",containment=True)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc5",containment=True)





# def plotEffectiveArea(zenChargeList):
# 	'''
# 	plots trigger efficiency in different zenith bins
# 	'''
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(nrows=1,ncols=1)
# 	ax = fig.add_subplot(gs[0])
# 	for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# 		zclist = [zc for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
# 		smtTrigEfficiencyList, globalTrigEfficiencyList = energyEfficiency(energyBins,zclist,"effArea")
# 		ax.plot(energyBins[:-1],smtTrigEfficiencyList,".",ls='-',lw = 2.5,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
# 	ax.set_xlabel(r"energy [GeV]", fontsize=24)
# 	ax.set_ylabel(r"effective area [$km^{2}$]", fontsize=24)
# 	# ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# 	ax.set_xscale('log')
# 	ax.grid(True,alpha=0.2)
# 	ax.legend(fontsize=12)
# 	plt.savefig(plotFolder+"/trigEffectiveArea.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# def rateZen(hdfFileList,zenBins,weighting,suffix):
# 	"""
# 	plots number of triggered events in different zenith angle bins
# 	"""
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	angleBins = []
# 	sta1_list = []
# 	sta3_list = []
# 	totalEvts_list = []
# 	sta1Frac_list = []
# 	sta3Frac_list = []
# 	for n,bin_start in enumerate(zenBins[:-1]):
# 		lowEdge = zenBins[n]
# 		highEdge = zenBins[n+1]
# 		totalEvts,sta1_evts,sta3_evts = getTriggeredEvents(hdfFileList,[lowEdge,highEdge],weighting)
# 		angleBins.append((lowEdge+highEdge)/2.0)
# 		sta1_list.append(sta1_evts)
# 		sta3_list.append(sta3_evts)
# 		totalEvts_list.append(totalEvts)
# 	ax.plot(angleBins,np.asarray(sta1_list)/trigWindow,"o-",c=next(colorsIter),label="STA1",alpha=1)
# 	ax.plot(angleBins,np.asarray(sta3_list)/trigWindow,"o-",c=next(colorsIter),label="STA3",alpha=1)
# 	# ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
# 	ax.set_ylabel(r"rate [Hz]", fontsize=20)
# 	ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
# 	ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
# 	# ax.set_yscale('log')
# 	ax.set_xlim(0,65)
# 	ax.set_title("total rate STA1 {:.2f} Hz STA3 {:.2f} Hz".format(np.sum(sta1_list)/trigWindow,np.sum(sta3_list)/trigWindow),fontsize=12)
# 	ax.grid(True,alpha=0.2)
# 	ax.legend(fontsize=14)
# 	# ax.legend(fontsize=14,ncol=2)
# 	plt.savefig(plotFolder+"/"+str(suffix)+"Wts"+str(weighting)+"triggerRate_zen.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()
# # rateZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,"HLCVEM")
# # rateZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,"HLCVEM")





# def delta_t_events_(hdfFile,events,key):
# 	'''gets delta_t out of given hdffile
# 	belonging to given events
# 	'''
# 	df_delta_t = pd.read_hdf(hdfFile,key=key)
# 	delta_t_selected = []
# 	events_t = df_delta_t["Event"].values
# 	delta_t = df_delta_t["item"].values
# 	for jevent,jdelta_t in zip(events_t,delta_t):
# 		if jevent in events:
# 			delta_t_selected.append(jdelta_t)
# 	return np.asarray(delta_t_selected)
# def delta_t_zenith_(hdfFile,zenLim,key):
# 	sel_evts = getEventsZenith_(hdfFile,zenLim)
# 	return delta_t_events_(hdfFile,sel_evts,key)

# def delta_t_zenith(hdfFileList,zenLim,key):
# 	delta_t_list = []
# 	for ihdf in hdfFileList:
# 		delta_t = delta_t_zenith_(ihdf,zenLim,key)
# 		delta_t_list = np.concatenate((delta_t_list,delta_t))
# 	return delta_t_list

# def delta_t_hist_zen_bins(hdfFileList,zenBins,suffix,key):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     delta_ts = getVectorItem(hdfFileList,key)
#     lastBin = min(3000,max(delta_ts))
#     delta_t_bins = np.linspace(min(delta_ts),lastBin,300)
#     for n,bin_start in enumerate(zenBins[:-1]):
#     	lowEdge = zenBins[n]
#     	highEdge = zenBins[n+1]
#     	delta_t_this = delta_t_zenith(hdfFileList,zenLim=[lowEdge,highEdge],key=key)
#     	print("plotting hist",lowEdge,highEdge,len(delta_t_this))
#     	ax.hist(delta_t_this,bins=delta_t_bins,
#     		histtype="step",label=r"$\theta$ = {0:.0f}$^{{\circ}}$-{1:.0f}$^{{\circ}}$, {2:d} evts".format(lowEdge,highEdge,len(delta_t_this)),color=next(colorsIter),lw=2.5,alpha=1)
#     ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
#     ax.set_xlabel(r"$\Delta_t$ [ns]", fontsize=20)
#     ax.set_ylabel(r"count", fontsize=20)
#     ax.set_yscale('log')
#     # ax.set_xscale('log')
#     # ax.set_xlim(None,2000)
#     ax.grid(True,alpha=0.2)
#     ax.set_title(key,fontsize=12)
#     ax.legend(fontsize=14)
#     # ax.legend(fontsize=14,ncol=2)
#     plt.savefig(plotFolder+"/"+str(suffix)+"delta_t_zenBins.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # delta_t_hist_zen_bins(hdf5NullList,[0,10,20,30,40,50,65.1],"HLCVEM","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t")
# # delta_t_hist_zen_bins(hdf5NullList,[0,10,20,30,40,50,65.1],"SLCVEM","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_delta_t")

# def getZenith_(hdfFile):
# 	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
# 	return np.rad2deg(df_MCPrimary["zenith"].values)
# def getZenith(hdfFileList):
# 	zenithList = np.array([])
# 	for ihdf in hdfFileList:
# 		izenith = getZenith_(ihdf)
# 		print("izenith")
# 		if len(izenith)>0:
# 			zenithList = np.concatenate((zenithList,izenith))
# 	return zenithList

# def getZenithBug_(hdfFile):
# 	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
# 	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
# 	return np.rad2deg(df_MCPrimary["zenith"].values)
# def getZenith(hdfFileList):
# 	zenithList = np.array([])
# 	for ihdf in hdfFileList:
# 		izenith = getZenith_(ihdf)
# 		print("izenith")
# 		if len(izenith)>0:
# 			zenithList = np.concatenate((zenithList,izenith))
# 	return zenithList

# def histZen(hdfFileList,suffix):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     zenithList = getZenith(hdfFileList)
#     zenithBins = np.linspace(min(zenithList),max(zenithList),100)
#     ax.hist(zenithList,bins=zenithBins,histtype="step",lw=2.5,alpha=1)
#     ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
#     ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
#     ax.set_ylabel(r"count", fontsize=20)
#     ax.set_yscale('log')
#     # ax.set_xlim(None,2000)
#     ax.grid(True,alpha=0.2)
#     # ax.set_title(fontsize=12)
#     ax.legend(fontsize=14)
#     plt.savefig(plotFolder+"/"+str(suffix)+"zenithDist.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # histZen(hdf5NullList,"hist")

# def applyWeight(hitList,weights):
# 	return np.multiply(hitList,weights)

# SLCHits = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalHit")
# # SLCHitsW = applyWeight(SLCHits,adjustedWeights)
# SLCDuration = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration")
# SLC_Q = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalCharge")
# # SLC_QW = applyWeight(SLC_Q,adjustedWeights)

# HLCHits = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalHit")
# # HLCHitsW = applyWeight(HLCHits,adjustedWeights)
# HLCDuration = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeHitTimeDuration")
# HLC_Q = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalCharge")
# # HLC_QW = applyWeight(HLC_Q,adjustedWeights)




# def getTimeCharge(key,hdf5List):
#     timeList = np.array([])
#     chargeList = np.array([])
#     for ihdf in hdf5List:
#         time,charge = getTimeCharge_(key,ihdf)
#         timeList = np.concatenate((timeList,time))
#         chargeList = np.concatenate((chargeList,charge))
#     return timeList, chargeList

# def getTimeCharge_(key,hdfFile):
#     '''takes HLC/SLC pulses key and hdf file, returns time and charge
#     '''
#     dataT = pd.read_hdf(hdfFile,key=key)
#     return dataT["time"].values,dataT["charge"].values


# def plotTimeHistogram(time,suffix,key):
#     time = np.log10(time)
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     hitBins = np.linspace(min(time),max(time),200)
#     ax.hist(time,bins=hitBins,histtype="step",lw=2.5,alpha=1)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#     ax.set_xlabel(r"log10(time [ns])", fontsize=20)
#     ax.set_ylabel(r"count", fontsize=20)
#     ax.set_yscale('log')
#     ax.grid(True,alpha=0.2)
#     ax.set_title(key,fontsize=16)
#     # ax.legend(fontsize=20)
#     plt.savefig(plotFolder+"/"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
#     plt.close()

# def plotTimeChargeScatter(time,charge,suffix,key):
#     time = np.log10(time)
#     charge = np.log10(charge)
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     ax.scatter(time,charge,s=10)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#     ax.set_xlabel(r"log10(time [ns])", fontsize=20)
#     ax.set_ylabel(r"log10(charge)", fontsize=20)
# #     ax.set_yscale('log')
#     ax.set_ylim(-4.5,4.5)
#     ax.set_title(key,fontsize=16)
#     ax.grid(True,alpha=0.2)
#     plt.savefig(plotFolder+"/scatterChargeTime"+str(suffix)+".png",transparent=False,bbox_inches='tight')
#     plt.close()
    


# def scatter_hist_(x, y, ax, ax_histx, ax_histy):
#     # no labels
#     ax_histx.tick_params(axis="x", labelbottom=False)
#     ax_histy.tick_params(axis="y", labelleft=False)

#     # the scatter plot:
#     ax.scatter(x, y)

#     # now determine nice limits by hand:
#     binwidth = 0.25
#     xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
#     lim = (int(xymax/binwidth) + 1) * binwidth

#     bins = np.arange(-lim, lim + binwidth, binwidth)
#     ax_histx.hist(x, bins=bins)
#     ax_histy.hist(y, bins=bins, orientation='horizontal')

# def scatter_hist(x,y,suffix):
#     fig = plt.figure(figsize=(8, 5))
#     gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7),
#                           left=0.1, right=0.9, bottom=0.1, top=0.9,
#                           wspace=0.05, hspace=0.05)
#     ax = fig.add_subplot(gs[1, 0])
#     ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
#     ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
#     # use the previously defined function
#     scatter_hist_(x, y, ax, ax_histx, ax_histy)
#     plt.savefig(plotFolder+"/scatterChargeTimeHist"+str(suffix)+".png",transparent=False,bbox_inches='tight')
# #     plt.show()
#     plt.close()

# def plotScatterHist(hdfFileList,key,suffix):
#     timeList,chargeList = getTimeCharge(key=key,hdf5List=hdfFileList)
#     plotTimeHistogram(timeList,suffix,key)
#     plotTimeChargeScatter(timeList,chargeList,suffix,key)
# #     scatter_hist(timeList,chargeList,suffix)
# #     print("checking length",len(timeList),len(chargeList),len(hdfFileList))


# # plotScatterHist(hdf5NullList,key="OfflineIceTopSLCVEMPulses",suffix="SLCVEMTime")
# # plotScatterHist(hdf5NullList,key="OfflineIceTopHLCVEMPulses",suffix="HLCVEMTime")
# # plotScatterHist(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge",suffix="SLCVEMTimeCleaned")
# # plotScatterHist(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge",suffix="HLCVEMTimeCleaned")
# # dataT = pd.read_hdf(hdf5NullList[0],key="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration")


# def getTankHitDuration(keyDuration,keyHit,hdf5List):
#     hitDurationList = np.array([])
#     hitList = np.array([])
#     for ihdf in hdf5List:
#         ihitDuration,ihit = getTankHitDuration_(keyDuration,keyHit,ihdf)
#         hitDurationList = np.concatenate((hitDurationList,ihitDuration))
#         hitList = np.concatenate((hitList,ihit))
#     return hitDurationList,hitList

# def getTankHitDuration_(keyDuration,keyHit,hdfFile):
#     '''takes HLC/SLC pulses key and hdf file, returns time and charge
#     '''
#     dataDuration = pd.read_hdf(hdfFile,key=keyDuration)
#     dataHit = pd.read_hdf(hdfFile,key=keyHit)
#     return dataDuration["value"].values,dataHit["value"].values


# def plotHitDuration(hit,HitDuration,prefix):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     ax.scatter(hit,HitDuration,s=10)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#     ax.set_xlabel(r"tank hit", fontsize=20)
#     ax.set_ylabel(r"hit duration [ns]", fontsize=20)
#     # ax.set_yscale('log')
#     ax.set_ylim(0.1,None)
#     ax.grid(True,alpha=0.2)
#     plt.savefig(plotFolder+"/scatterHitDuration"+str(prefix)+".png",transparent=False,bbox_inches='tight')
#     plt.close()




# # slcHitDuration,slcHit = getTankHitDuration(keyDuration="OfflineIceTopSLCTankPulsesHitTimeDuration",keyHit="OfflineIceTopSLCTankPulsesTotalHit",hdf5List=hdf5NullList)
# # hlcHitDuration,hlcHit = getTankHitDuration(keyDuration="OfflineIceTopHLCTankPulsesHitTimeDuration",keyHit="OfflineIceTopHLCTankPulsesTotalHit",hdf5List=hdf5NullList)
# # plotHitDuration(slcHit,slcHitDuration,prefix="SLC")
# # plotHitDuration(hlcHit,hlcHitDuration,prefix="HLC")
# # slcHitDuration,slcHit = getTankHitDuration(keyDuration="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration",keyHit="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalHit",hdf5List=hdf5NullList)
# # hlcHitDuration,hlcHit = getTankHitDuration(keyDuration="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeHitTimeDuration",keyHit="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalHit",hdf5List=hdf5NullList)
# # plotHitDuration(slcHit,slcHitDuration,prefix="SLCClean")
# # plotHitDuration(hlcHit,hlcHitDuration,prefix="HLCClean")

# print("slc hit w",SLCHits)
# # print("slc hit w",SLCHitsW)
# # plotHitDuration(SLCHitsW,SLCDuration,"SLCWeightHit")


# def hitsWZenith_(hdfFile,key):
# 	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
# 	zenithList = np.rad2deg(df_MCPrimary["zenith"].values)
# 	hitList = pd.read_hdf(hdfFile,key=key)["value"].values
# 	return hitList,zenithList

# def hitsWZenith(hdfFileList,key):
# 	hitList = np.array([])
# 	zenithList = np.array([])
# 	for ihdf in hdfFileList:
# 		hits,zeniths = hitsWZenith_(ihdf,key)
# 		hitList = np.concatenate((hitList,hits))
# 		zenithList = np.concatenate((zenithList,zeniths))
# 	return hitList,zenithList
# def plotScatterHitsZenith(hdfFileList,key,prefix):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     hitList,zenithList = hitsWZenith(hdfFileList,key)
#     ax.scatter(zenithList,hitList,s=10)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#     ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
#     ax.set_ylabel(r"" + str(prefix)+" hits [ns]", fontsize=20)
#     # ax.set_yscale('log')
#     ax.set_title(key,fontsize=16)
#     ax.set_ylim(0.1,None)
#     ax.grid(True,alpha=0.2)
#     plt.savefig(plotFolder+"/scatterHitzenith"+str(prefix).replace(" ","")+".png",transparent=False,bbox_inches='tight')
#     plt.close()

# # plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopHLCVEMPulsesTotalHit",prefix="HLC VEM")
# # plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopSLCVEMPulsesTotalHit",prefix="SLC VEM")

# # plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalHit",prefix="Clean HLC VEM")
# # plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalHit",prefix="Clean SLC VEM")

# def eventsZenith(zenLim):
# 	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
# 	zeniths = dataT["zenith"].values


# def compareTwoWeights(evtList):
# 	frameWeight = [ievt.H4aWeight for ievt in evtList]
# 	calcWeight = [ievt.directWeight for ievt in evtList]
# 	print("length of two weights",len(frameWeight),len(calcWeight))
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	# hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
# 	hitBins = np.linspace(min(min(calcWeight),min(calcWeight)),max(max(calcWeight),max(calcWeight)),100)
# 	ax.hist(calcWeight,bins=hitBins,histtype="step",lw=2.5,label=r"Weight(calc)",alpha=1)
# 	ax.hist(frameWeight,bins=hitBins,histtype="step",lw=2.5,label="Weight(frame)",alpha=1)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
# 	ax.set_xlabel(r"weight", fontsize=20)
# 	ax.set_ylabel(r"count", fontsize=20)
# 	ax.set_yscale("log")
# 	# ax.set_ylim(None,10**5)
# 	# ax.set_xscale('log')
# 	ax.grid(True,alpha=0.2)
# 	# ax.set_title(key,fontsize=16)
# 	ax.legend(fontsize=10,ncol=3,loc="lower center")
# 	plt.savefig(plotFolder+"/weightCompare.pdf",transparent=False,bbox_inches='tight')
# 	plt.close()
# # compareTwoWeights(evtList)

# def stationRate(hdfFileList,key,prefix,hitBins,xlabel):
# 	stationHits = getValue(hdfFileList,key)
# 	weights = getValue(hdfFileList,"H4aWeight")
# 	triggerSTA3 = getValue(hdfFileList,"ITSMTTriggered")
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	# hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
# 	ax.hist(stationHits,bins=hitBins,histtype="step",weights=[w*t for w,t in zip(weights,triggerSTA3)],lw=2.5,label=r"stations",alpha=1)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=12)
# 	ax.tick_params(axis='both',which='major',length=5)
# 	ax.tick_params(axis='both',which='minor',length=3)
# 	ax.set_xlabel(xlabel, fontsize=12)
# 	ax.set_ylabel(r"rate [Hz]", fontsize=12)
# 	ax.set_yscale("log")
# 	ticks = hitBins
# 	ax.set_xticks(ticks[::10])
# 	ax.set_xticks(ticks[::2],True)
# 	ax.set_xlim(0,None)
# 	# ax.set_xscale('log')
# 	# ax.grid(True,alpha=0.2)
# 	# ax.set_title(key,fontsize=16)
# 	# ax.legend(fontsize=10,ncol=3,loc="lower center")
# 	plt.savefig(plotFolder+"/Rate"+str(prefix)+".pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # stationRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_hitStations",prefix="HLCTankStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (HLC Tank)")
# # stationRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_hitStations",prefix="SLCTankStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (SLC Tank)")
# # stationRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_hitTanks",prefix="HLCTankTank",hitBins=np.linspace(0,160,161),xlabel=r"$N_{tanks}$ per event (HLC Tank)")
# # stationRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_hitTanks",prefix="SLCTankTank",hitBins=np.linspace(0,80,81),xlabel=r"$N_{tanks}$ per event (SLC Tank)")

# # stationRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_hitStations",prefix="HLCVEMStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (HLC VEM)")
# # stationRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_hitStations",prefix="SLCVEMStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (SLC VEM)")
# # stationRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_hitTanks",prefix="HLCVEMTank",hitBins=np.linspace(0,160,161),xlabel=r"$N_{tanks}$ per event (HLC VEM)")
# # stationRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_hitTanks",prefix="SLCVEMTank",hitBins=np.linspace(0,80,81),xlabel=r"$N_{tanks}$ per event (SLC VEM)")




# def chargeRate(hdfFileList,key,prefix,hitBins,xlabel):
# 	Qtot = np.log10(getValue(hdfFileList,key))
# 	weights = getValue(hdfFileList,"H4aWeight")
# 	triggerSTA3 = getValue(hdfFileList,"ITSMTTriggered")
# 	fig = plt.figure(figsize=(8,5))
# 	gs = gridspec.GridSpec(ncols=1,nrows=1)
# 	ax = fig.add_subplot(gs[0])
# 	# hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
# 	ax.hist(Qtot,bins=hitBins,histtype="step",weights=[w*t for w,t in zip(weights,triggerSTA3)],lw=2.5,label=r"stations",alpha=1)
# 	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
# 	ax.tick_params(axis='both',which='both', direction='in', labelsize=12)
# 	ax.tick_params(axis='both',which='major',length=5)
# 	ax.tick_params(axis='both',which='minor',length=3)
# 	ax.set_xlabel(xlabel, fontsize=12)
# 	ax.set_ylabel(r"rate [Hz]", fontsize=12)
# 	ax.set_yscale("log")
# 	ticks = hitBins
# 	ax.set_xticks(ticks[::10])
# 	ax.set_xticks(ticks[::2],True)
# 	# ax.set_xlim(0,None)
# 	# ax.set_xscale('log')
# 	# ax.grid(True,alpha=0.2)
# 	# ax.set_title(key,fontsize=16)
# 	# ax.legend(fontsize=10,ncol=3,loc="lower center")
# 	plt.savefig(plotFolder+"/chargeRate"+str(prefix)+".pdf",transparent=False,bbox_inches='tight')
# 	plt.close()

# # chargeRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalCharge",prefix="HLCTankCharge",hitBins=np.linspace(0,5,81),xlabel=r"log($Q_{tot}/VEM$) per event (HLC tank)")
# # chargeRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalCharge",prefix="SLCTankCharge",hitBins=np.linspace(-2,4,81),xlabel=r"log($Q_{tot}/VEM$) per event (SLC tank)")
# # chargeRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalCharge",prefix="HLCVEMCharge",hitBins=np.linspace(0,5,81),xlabel=r"log($Q_{tot}/VEM$) per event (HLC VEM)")
# # chargeRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalCharge",prefix="SLCVEMCharge",hitBins=np.linspace(-2,4,81),xlabel=r"log($Q_{tot}/VEM$) per event (SLC VEM)")




def SLCRate(evtList):
	# print([ievt.zenith for ievt in evtList if 0<np.rad2deg(ievt.zenith)<60][:100])
	evtList = [ievt for ievt in evtList if 60<np.rad2deg(ievt.zenith)<65]
	nSLCList = [ievt.nSLC for ievt in evtList]
	nSLCList_sta3 = [ievt.nSLC for ievt in evtList if abs(ievt.ITSMTTriggered-1)<0.01]
	nSLCList_notsta3 = [ievt.nSLC for ievt in evtList if abs(ievt.ITSMTTriggered-1)>0.01]
	weights = [ievt.H4aWeight for ievt in evtList]
	weights_sta3 = [ievt.H4aWeight for ievt in evtList if abs(ievt.ITSMTTriggered-1)<0.01]
	weights_notsta3 = [ievt.H4aWeight for ievt in evtList if abs(ievt.ITSMTTriggered-1)>0.01]
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(ncols=1,nrows=1)
	ax = fig.add_subplot(gs[0])
	hitBins = np.linspace(0,40,41)
	n, bins, patches=ax.hist(nSLCList,bins=hitBins,histtype="step",weights=weights,lw=2.5,label=r"nSLC",alpha=1)
	n1, bins1, patches1=ax.hist(nSLCList_sta3,bins=hitBins,histtype="step",weights=weights_sta3,lw=2.5,label=r"nSLC sta3",alpha=1)
	n2, bins2, patches2=ax.hist(nSLCList_notsta3,bins=hitBins,histtype="step",weights=weights_notsta3,lw=2.5,label=r"nSLC not trig",alpha=1)
	# print(n,bins)
	# print(n1,bins1)
	# print(n2,bins2)
	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=12)
	ax.tick_params(axis='both',which='major',length=5)
	ax.tick_params(axis='both',which='minor',length=3)
	ax.set_xlabel("n_SLC", fontsize=12)
	ax.set_ylabel(r"rate [Hz]", fontsize=12)
	ax.set_yscale("log")
	ticks = hitBins
	ax.set_xticks(ticks[::10])
	ax.set_xticks(ticks[::1],True)
	# ax.set_xlim(0,None)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.4)
	# ax.set_title(key,fontsize=16)
	ax.legend(fontsize=10,ncol=2,loc="lower center")
	plt.savefig(plotFolder+"/SLCRate.pdf",transparent=False,bbox_inches='tight')
	plt.close()

# SLCRate(evtList)




