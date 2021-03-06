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
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap

import numpy as np

from weighting import GetWeight

basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTest/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTest/"
basePathOfficial = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanOfficialL2/"
# basePathOfficial = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanOfficial/"
hdf5NullListP = sorted(glob.glob(basePath+"p*Clean*.hdf5"))
hdf5NullListHe = sorted(glob.glob(basePath+"He*Clean*.hdf5"))
hdf5NullListO = sorted(glob.glob(basePath+"O*Clean*.hdf5"))
hdf5NullListFe = sorted(glob.glob(basePath+"Fe*Clean*.hdf5"))
hdf5NullList = np.concatenate((hdf5NullListP,hdf5NullListHe,hdf5NullListO,hdf5NullListFe))


###############################################################################
############################for l3
###################################

# hdf5NullListPOfficial = sorted(glob.glob(basePathOfficial+"2012_12360*.hdf5"))
# hdf5NullListHeOfficial = sorted(glob.glob(basePathOfficial+"2012_12630*.hdf5"))
# hdf5NullListOOfficial = sorted(glob.glob(basePathOfficial+"2012_12631*.hdf5"))
# hdf5NullListFeOfficial = sorted(glob.glob(basePathOfficial+"2012_12362*.hdf5"))

####################################################################################
##for L2
####################################################################################
hdf5NullListPOfficial = sorted(glob.glob(basePathOfficial+"010410*.hdf5"))
hdf5NullListHeOfficial = sorted(glob.glob(basePathOfficial+"011663*.hdf5"))
hdf5NullListOOfficial = sorted(glob.glob(basePathOfficial+"012605*.hdf5"))
hdf5NullListFeOfficial = sorted(glob.glob(basePathOfficial+"012362*.hdf5"))

hdf5NullListOfficial = np.concatenate((hdf5NullListPOfficial,hdf5NullListHeOfficial,hdf5NullListOOfficial,hdf5NullListFeOfficial))
print("official",len(hdf5NullListOfficial),hdf5NullListOfficial[:5])

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"


colorsIter = iter(['#1f77b4','#9467bd','#e377c2', '#ff7f0e','#bcbd22', '#2ca02c', '#8c564b', '#e377c2', '#7f7f7f',  '#17becf','#d62728',"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray'])

def nCorFiles(hdfFileList):
	nfiles = 0
	for ihdf in hdfFileList:
		dataT = pd.read_hdf(ihdf,key="I3EventHeader")
		nfiles += len(set(dataT["Run"].values))
	return nfiles

nfilesP = nCorFiles(hdf5NullListP)
nfilesHe = nCorFiles(hdf5NullListHe)
nfilesO = nCorFiles(hdf5NullListO)
nfilesFe = nCorFiles(hdf5NullListFe)
print("no of files",nfilesP,nfilesHe,nfilesO,nfilesFe)
nfilesP = nCorFiles(hdf5NullListPOfficial)
nfilesHe = nCorFiles(hdf5NullListHeOfficial)
nfilesO = nCorFiles(hdf5NullListOOfficial)
nfilesFe = nCorFiles(hdf5NullListFeOfficial)
print("no of files",nfilesP,nfilesHe,nfilesO,nfilesFe)


trigWindow = 10**(-6) # in ns

def getValue_(hdfFile,key):
	'''extract value of given key in hdf file 
	as single array
	'''
	dataT = pd.read_hdf(hdfFile,key=key)
	return dataT["value"].values

def getValueFilt_(hdfFile,key,filtKey):
	'''extract value of given key in hdf file 
	as single array and filter
	'''
	keyValues = getValue_(hdfFile,key)
	filterValues = getValue_(hdfFile,filtKey)
	# print("zeros and ones",len(filterValues),np.count_nonzero(filterValues==0),np.count_nonzero(filterValues))
	return keyValues[abs(np.asarray(filterValues-1)) < 0.001]

def getValueFilt(hdfFileList,key,filtKey):
	hitList = np.array([])
	for ihdf in hdfFileList:
		ihit = getValueFilt_(ihdf,key,filtKey)
		hitList = np.concatenate((hitList,ihit))
	return hitList

def getValue(hdfFileList,key):
	hitList = np.array([])
	for ihdf in hdfFileList:
		ihit = getValue_(ihdf,key)
		hitList = np.concatenate((hitList,ihit))
	return hitList

def getVectorItem(hdfFileList,key):
	hitDuration = np.array([])
	for ihdf in hdfFileList:
		dt_hdf = pd.read_hdf(ihdf,key=key)["item"].values
		hitDuration = np.concatenate((hitDuration,dt_hdf))
	return hitDuration

def getZenith_(hdfFile):
	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
	return dataT["zenith"].values

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

def getZenithTypeEnergyFilt_(hdfFile,filtKey):
	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
	# filtBool = getValue_(hdfFile,filtKey)
	# filtBool = getValue_(hdfFile,filtKey) == 1
	filtBool = abs(np.asarray(getValue_(hdfFile,filtKey) - 1.0)) < 0.001
	return dataT["zenith"].values[filtBool], dataT["type"].values[filtBool], dataT["energy"].values[filtBool]

def getZenithTypeEnergyFilt(hdfFileList,filtKey):
	ptypeList = np.array([])
	energyList = np.array([])
	zenList = np.array([])
	for ihdf in hdfFileList:
		zen,ptype,energy = getZenithTypeEnergyFilt_(ihdf,filtKey)
		zenList = np.concatenate((zenList,zen))
		ptypeList = np.concatenate((ptypeList,ptype))
		energyList = np.concatenate((energyList,energy))
	return zenList,ptypeList,energyList


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


# zenith,ptype,energy = getZenithTypeEnergy(hdf5NullList)




# weights = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zenith,energy,ptype)
# weithtsOne = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zenith[-1],energy[-1],ptype[-1])
# print("weightCompare",weights[-1],weithtsOne,len(weights))
# adjustedWeights = []

# for iweight,ienergy in zip(weights,energy):
# 	if ienergy >= 10**7:
# 		adjustedWeights.append(iweight*10)
# 	elif ienergy < 10**6:
# 		adjustedWeights.append(iweight*10)
# 	else:
# 		adjustedWeights.append(iweight)


def plotEnergySpectraCompare(hdfFileList,hdfFileListOfficial,filtKey,filtKeyOfficial,yscale,suffix,energyScale):
	zen,ptype,energy = getZenithTypeEnergy(hdfFileList)
	print("entries before filt",len(energy))
	zen,ptype,energy = getZenithTypeEnergyFilt(hdfFileList,filtKey)
	print("entries after filt",len(energy))
	weights = getValue(hdfFileList,"H4aWeight")
	print("weight entries before filt",len(weights))
	weights = getValueFilt(hdfFileList,"H4aWeight",filtKey)
	print("weight entries after filt",len(weights))
	zenOfficial,ptypeOfficial,energyOfficial = getZenithTypeEnergyFilt(hdfFileListOfficial,filtKeyOfficial)
	# zenOfficial,ptypeOfficial,energyOfficial = getZenithTypeEnergyFilt(hdfFileListOfficial,"SDST_IceTopSTA3_12_filter")
	print("entries after filt",len(energyOfficial))
	weightsOfficial = getValueFilt(hdfFileListOfficial,"H4aWeight",filtKey=filtKeyOfficial)
	# weightsOfficial = getValueFilt(hdfFileListOfficial,"H4aWeight",filtKey="SDST_IceTopSTA3_12_filter")
	print("weight entries after filt",len(weightsOfficial))

	print("length check",len(weightsOfficial),len(energyOfficial))
	print("length check",weightsOfficial[:5],energyOfficial[:5])
	adjustedWeights = []
	for iweight,ienergy in zip(weights,energy):
		if ienergy >= 10**7:
			# adjustedWeights.append(iweight*0.78)
			adjustedWeights.append(iweight)
		elif ienergy < 10**6:
			adjustedWeights.append(iweight)
		else:
			adjustedWeights.append(iweight)
	weights = adjustedWeights
	energy = np.log10(energy)
	energyOfficial = np.log10(energyOfficial)
	# print("lengths",len(weights),weights,len(energy))
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(ncols=1,nrows=1)
	ax = fig.add_subplot(gs[0])
	hitBins = np.linspace(5.0,8.0,31)
	hist,binEdge = np.histogram(energy,hitBins,weights=weights)
	histOfficial,binEdge = np.histogram(energyOfficial,hitBins,weights=weightsOfficial)
	print("sum of official rate",sum(histOfficial))
	print("sum of simulated rate",sum(hist))
	ax.text(0.05, 0.95, r"{} official rate:{:.1f} Hz, sim. rate:{:.1f} Hz".format(filtKeyOfficial,sum(histOfficial),sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
	binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
	H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
	HOfficial = [h*(10**E)**energyScale for h,E in zip(histOfficial,binCenter)]
	######################################################################################
	# ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
	# ax.hist(energyOfficial,bins=hitBins,histtype="step",lw=2.5,label=r"after weighting",alpha=1)
	####################################################################################
	ax.step(binCenter,H,"-",where="mid",lw=2.5,label="my simulation",alpha=0.6)
	ax.step(binCenter,HOfficial,"-",where="mid",lw=2.5,label="official simulation",alpha=0.6)
	####################################################################################
	# ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
	ax.set_xlabel(r"log10(energy [GeV])", fontsize=20)
	ax.set_ylabel(r"$E^{0:.1f}$ rate [Hz]".format(energyScale), fontsize=20)
	ax.set_yscale(yscale)
	# ax.set_ylim(10**10,10**11.5)
	ax.set_ylim(10**6,10**12)
	# ax.set_xlim(6.5,7.5)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.set_title(key,fontsize=16)
	ax.legend(fontsize=10,ncol=3,loc="lower center")
	plt.savefig(plotFolder+"/energySpec"+str(suffix)+"Compare.pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotEnergySpectraCompare(hdf5NullList,hdf5NullListOfficial,"IceTopSTA5_13_filter","IceTopSTA5_12_filter","log","fluxLogSTA5",energyScale=1.8)
plotEnergySpectraCompare(hdf5NullList,hdf5NullListOfficial,"SDST_IceTopSTA3_13_filter","SDST_IceTopSTA3_12_filter","log","fluxLogSTA3",energyScale=1.8)



def plotEnergyFlux(hdfFileList,yscale,suffix):
	zen,ptype,energy = getZenithTypeEnergy(hdfFileList)
	weights = getValue(hdfFileList,"H4aWeight")
	adjustedWeights = []
	for iweight,ienergy in zip(weights,energy):
		if ienergy >= 10**7:
			adjustedWeights.append(iweight/1.0)
		elif ienergy < 10**6:
			adjustedWeights.append(iweight)
		else:
			adjustedWeights.append(iweight)
	weights = adjustedWeights
	triggerSTA3 = getValue(hdfFileList,"ITSMTTriggered")
	energy = np.log10(energy)
	# print("lengths",len(weights),weights,len(energy))
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(ncols=1,nrows=1)
	ax = fig.add_subplot(gs[0])
	hitBins = np.linspace(5.0,8.0,31)
	hist,binEdge = np.histogram(energy,hitBins,weights=[w*t for w,t in zip(weights,triggerSTA3)])
	binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
	H = [h*(10**E)**0 for h,E in zip(hist,binCenter)]
	# ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
	ax.step(binCenter,H,"-",where="mid",lw=2.5,label="after weighting",alpha=1)
	# ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
	ax.set_xlabel(r"log10(energy [GeV])", fontsize=20)
	ax.set_ylabel(r"$E^{1.8}$ rate [Hz]", fontsize=20)
	ax.set_yscale(yscale)
	# ax.set_ylim(None,10**5)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.set_title(key,fontsize=16)
	ax.legend(fontsize=10,ncol=3,loc="lower center")
	plt.savefig(plotFolder+"/energySpec"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()
# plotEnergyFlux(energy,weights,adjustedWeights,"log","flux")
plotEnergyFlux(hdf5NullList,"linear","fluxlinear")
plotEnergyFlux(hdf5NullList,"log","fluxLog")



HLCStationDuration = getVectorItem(hdf5NullList,"OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t")
def delta_t_hist(duration,suffix,key):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    hitBins = np.linspace(min(duration),max(duration),200)
    ax.hist(duration,bins=hitBins,histtype="step",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"$\Delta_t$ [ns]", fontsize=20)
    ax.set_ylabel(r"count", fontsize=20)
    ax.set_yscale('log')
    ax.set_xlim(None,2000)
    ax.grid(True,alpha=0.2)
    ax.set_title(key,fontsize=12)
    # ax.legend(fontsize=20)
    plt.savefig(plotFolder+"/"+str(suffix)+"delta_t.pdf",transparent=False,bbox_inches='tight')
    plt.close()
# delta_t_hist(HLCStationDuration,"HLCVEM","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t")

def getEventsZenith_(hdfFile,zenLim):
	'''get events belonging to given zenith bin:
	zenith angles in degree
	'''	
	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
	events_MC = df_MCPrimary["Event"].values
	zeniths_MC = df_MCPrimary["zenith"].values
	selEvents = []
	for ievent,izenith in zip(events_MC,zeniths_MC):
		if np.rad2deg(izenith) >= zenLim[0] and np.rad2deg(izenith) < zenLim[1]:
			selEvents.append(ievent)
	return list(set(selEvents))

def getTriggeredEvents_(hdfFile,zenLim,weighting):
	"""
	returns total events sta1 and sta3 trigered events
	"""
	STA1Trigger_df = pd.read_hdf(hdfFile,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_isSTA1")	
	STA3Trigger_df = pd.read_hdf(hdfFile,key="ITSMTTriggered")
	selEvents = getEventsZenith_(hdfFile,zenLim)
	eventList = STA1Trigger_df["Event"].values
	totalEvts = np.ones(len(eventList))
	if weighting == True:
		# weights,_ = weightCalc(hdfFile)
		weights = getValue_(hdfFile,"H4aWeight")
		totalEvts = np.multiply(totalEvts,weights)
		STA1Triggers = np.multiply(STA1Trigger_df["value"].values,weights)
		STA3Triggers = np.multiply(STA3Trigger_df["value"].values,weights)
	else:
		STA1Triggers = STA1Trigger_df["value"].values
		STA3Triggers = STA3Trigger_df["value"].values
	STA1Triggers_select = []
	STA3Triggers_select = []
	totalEvts_select = []
	for ievt,ista1,ista3,isEvt in zip(eventList,STA1Triggers,STA3Triggers,totalEvts):
		if ievt in selEvents:
			STA1Triggers_select.append(ista1)
			STA3Triggers_select.append(ista3)
			totalEvts_select.append(isEvt)
	return np.sum(totalEvts),np.sum(STA1Triggers_select),np.sum(STA3Triggers_select)
	
def getTriggeredEvents(hdfFileList,zenLim,weighting):
	"""
	sums total sta1 and sta3 events from all hdfFileList files
	"""
	totalEvts = 0
	sta1_evts = 0
	sta3_evts = 0
	for ihdf in hdfFileList:
		n_total,n_sta1,n_sta3 = getTriggeredEvents_(ihdf,zenLim,weighting)
		totalEvts += n_total
		sta1_evts += n_sta1
		sta3_evts += n_sta3
	return totalEvts,sta1_evts,sta3_evts

def trigZen(hdfFileList,zenBins,fraction,weighting,suffix):
	"""
	plots number of triggered events in different zenith angle bins
	"""
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(ncols=1,nrows=1)
	ax = fig.add_subplot(gs[0])
	angleBins = []
	sta1_list = []
	sta3_list = []
	totalEvts_list = []
	sta1Frac_list = []
	sta3Frac_list = []
	for n,bin_start in enumerate(zenBins[:-1]):
		lowEdge = zenBins[n]
		highEdge = zenBins[n+1]
		totalEvts,sta1_evts,sta3_evts = getTriggeredEvents(hdfFileList,[lowEdge,highEdge],weighting)
		angleBins.append((lowEdge+highEdge)/2.0)
		sta1_list.append(sta1_evts)
		sta3_list.append(sta3_evts)
		totalEvts_list.append(totalEvts)
	if fraction == True:
		sta1Frac_list = np.asarray(sta1_list)/np.asarray(totalEvts_list)
		sta3Frac_list = np.asarray(sta3_list)/np.asarray(totalEvts_list)
		ax.plot(angleBins,sta1Frac_list,"o-",c=next(colorsIter),label="STA1",alpha=1)
		ax.plot(angleBins,sta3Frac_list,"o-",c=next(colorsIter),label="STA3",alpha=1)
		ax.set_ylabel(r"fraction", fontsize=20)
	elif fraction == False:
	    ax.plot(angleBins,sta1_list,"o-",c=next(colorsIter),label="STA1",alpha=1)
	    ax.plot(angleBins,sta3_list,"o-",c=next(colorsIter),label="STA3",alpha=1)
	    # ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
	    ax.set_ylabel(r"rate [Hz]", fontsize=20)
	ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
	ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
	# ax.set_yscale('log')
	ax.set_title("total rate STA1 {:.2f} Hz STA3 {:.2f} Hz".format(np.sum(sta1_list),np.sum(sta3_list)),fontsize=12)
	ax.set_xlim(0,65)
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=14)
	# ax.legend(fontsize=14,ncol=2)
	plt.savefig(plotFolder+"/"+str(suffix)+"Frac"+str(fraction)+"Wts"+str(weighting)+"trigger_zen.pdf",transparent=False,bbox_inches='tight')
	plt.close()

# trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,False,"HLCVEM")
# trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],False,False,"HLCVEM")
# trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,True,"HLCVEM")
# trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],False,True,"HLCVEM")

def rateZen(hdfFileList,zenBins,weighting,suffix):
	"""
	plots number of triggered events in different zenith angle bins
	"""
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(ncols=1,nrows=1)
	ax = fig.add_subplot(gs[0])
	angleBins = []
	sta1_list = []
	sta3_list = []
	totalEvts_list = []
	sta1Frac_list = []
	sta3Frac_list = []
	for n,bin_start in enumerate(zenBins[:-1]):
		lowEdge = zenBins[n]
		highEdge = zenBins[n+1]
		totalEvts,sta1_evts,sta3_evts = getTriggeredEvents(hdfFileList,[lowEdge,highEdge],weighting)
		angleBins.append((lowEdge+highEdge)/2.0)
		sta1_list.append(sta1_evts)
		sta3_list.append(sta3_evts)
		totalEvts_list.append(totalEvts)
	ax.plot(angleBins,np.asarray(sta1_list)/trigWindow,"o-",c=next(colorsIter),label="STA1",alpha=1)
	ax.plot(angleBins,np.asarray(sta3_list)/trigWindow,"o-",c=next(colorsIter),label="STA3",alpha=1)
	# ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
	ax.set_ylabel(r"rate [Hz]", fontsize=20)
	ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
	ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
	# ax.set_yscale('log')
	ax.set_xlim(0,65)
	ax.set_title("total rate STA1 {:.2f} Hz STA3 {:.2f} Hz".format(np.sum(sta1_list)/trigWindow,np.sum(sta3_list)/trigWindow),fontsize=12)
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=14)
	# ax.legend(fontsize=14,ncol=2)
	plt.savefig(plotFolder+"/"+str(suffix)+"Wts"+str(weighting)+"triggerRate_zen.pdf",transparent=False,bbox_inches='tight')
	plt.close()
# rateZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,"HLCVEM")
# rateZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,"HLCVEM")





def delta_t_events_(hdfFile,events,key):
	'''gets delta_t out of given hdffile
	belonging to given events
	'''
	df_delta_t = pd.read_hdf(hdfFile,key=key)
	delta_t_selected = []
	events_t = df_delta_t["Event"].values
	delta_t = df_delta_t["item"].values
	for jevent,jdelta_t in zip(events_t,delta_t):
		if jevent in events:
			delta_t_selected.append(jdelta_t)
	return np.asarray(delta_t_selected)
def delta_t_zenith_(hdfFile,zenLim,key):
	sel_evts = getEventsZenith_(hdfFile,zenLim)
	return delta_t_events_(hdfFile,sel_evts,key)

def delta_t_zenith(hdfFileList,zenLim,key):
	delta_t_list = []
	for ihdf in hdfFileList:
		delta_t = delta_t_zenith_(ihdf,zenLim,key)
		delta_t_list = np.concatenate((delta_t_list,delta_t))
	return delta_t_list

def delta_t_hist_zen_bins(hdfFileList,zenBins,suffix,key):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    delta_ts = getVectorItem(hdfFileList,key)
    lastBin = min(3000,max(delta_ts))
    delta_t_bins = np.linspace(min(delta_ts),lastBin,300)
    for n,bin_start in enumerate(zenBins[:-1]):
    	lowEdge = zenBins[n]
    	highEdge = zenBins[n+1]
    	delta_t_this = delta_t_zenith(hdfFileList,zenLim=[lowEdge,highEdge],key=key)
    	print("plotting hist",lowEdge,highEdge,len(delta_t_this))
    	ax.hist(delta_t_this,bins=delta_t_bins,
    		histtype="step",label=r"$\theta$ = {0:.0f}$^{{\circ}}$-{1:.0f}$^{{\circ}}$, {2:d} evts".format(lowEdge,highEdge,len(delta_t_this)),color=next(colorsIter),lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
    ax.set_xlabel(r"$\Delta_t$ [ns]", fontsize=20)
    ax.set_ylabel(r"count", fontsize=20)
    ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.set_xlim(None,2000)
    ax.grid(True,alpha=0.2)
    ax.set_title(key,fontsize=12)
    ax.legend(fontsize=14)
    # ax.legend(fontsize=14,ncol=2)
    plt.savefig(plotFolder+"/"+str(suffix)+"delta_t_zenBins.pdf",transparent=False,bbox_inches='tight')
    plt.close()
# delta_t_hist_zen_bins(hdf5NullList,[0,10,20,30,40,50,65.1],"HLCVEM","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t")
# delta_t_hist_zen_bins(hdf5NullList,[0,10,20,30,40,50,65.1],"SLCVEM","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_delta_t")

def getZenith_(hdfFile):
	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
	return np.rad2deg(df_MCPrimary["zenith"].values)
def getZenith(hdfFileList):
	zenithList = np.array([])
	for ihdf in hdfFileList:
		izenith = getZenith_(ihdf)
		print("izenith")
		if len(izenith)>0:
			zenithList = np.concatenate((zenithList,izenith))
	return zenithList

def getZenithBug_(hdfFile):
	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
	return np.rad2deg(df_MCPrimary["zenith"].values)
def getZenith(hdfFileList):
	zenithList = np.array([])
	for ihdf in hdfFileList:
		izenith = getZenith_(ihdf)
		print("izenith")
		if len(izenith)>0:
			zenithList = np.concatenate((zenithList,izenith))
	return zenithList

def histZen(hdfFileList,suffix):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    zenithList = getZenith(hdfFileList)
    zenithBins = np.linspace(min(zenithList),max(zenithList),100)
    ax.hist(zenithList,bins=zenithBins,histtype="step",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
    ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
    ax.set_ylabel(r"count", fontsize=20)
    ax.set_yscale('log')
    # ax.set_xlim(None,2000)
    ax.grid(True,alpha=0.2)
    # ax.set_title(fontsize=12)
    ax.legend(fontsize=14)
    plt.savefig(plotFolder+"/"+str(suffix)+"zenithDist.pdf",transparent=False,bbox_inches='tight')
    plt.close()
# histZen(hdf5NullList,"hist")

def applyWeight(hitList,weights):
	return np.multiply(hitList,weights)

SLCHits = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalHit")
# SLCHitsW = applyWeight(SLCHits,adjustedWeights)
SLCDuration = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration")
SLC_Q = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalCharge")
# SLC_QW = applyWeight(SLC_Q,adjustedWeights)

HLCHits = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalHit")
# HLCHitsW = applyWeight(HLCHits,adjustedWeights)
HLCDuration = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeHitTimeDuration")
HLC_Q = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalCharge")
# HLC_QW = applyWeight(HLC_Q,adjustedWeights)




def getTimeCharge(key,hdf5List):
    timeList = np.array([])
    chargeList = np.array([])
    for ihdf in hdf5List:
        time,charge = getTimeCharge_(key,ihdf)
        timeList = np.concatenate((timeList,time))
        chargeList = np.concatenate((chargeList,charge))
    return timeList, chargeList

def getTimeCharge_(key,hdfFile):
    '''takes HLC/SLC pulses key and hdf file, returns time and charge
    '''
    dataT = pd.read_hdf(hdfFile,key=key)
    return dataT["time"].values,dataT["charge"].values


def plotTimeHistogram(time,suffix,key):
    time = np.log10(time)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    hitBins = np.linspace(min(time),max(time),200)
    ax.hist(time,bins=hitBins,histtype="step",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"log10(time [ns])", fontsize=20)
    ax.set_ylabel(r"count", fontsize=20)
    ax.set_yscale('log')
    ax.grid(True,alpha=0.2)
    ax.set_title(key,fontsize=16)
    # ax.legend(fontsize=20)
    plt.savefig(plotFolder+"/"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plotTimeChargeScatter(time,charge,suffix,key):
    time = np.log10(time)
    charge = np.log10(charge)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    ax.scatter(time,charge,s=10)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"log10(time [ns])", fontsize=20)
    ax.set_ylabel(r"log10(charge)", fontsize=20)
#     ax.set_yscale('log')
    ax.set_ylim(-4.5,4.5)
    ax.set_title(key,fontsize=16)
    ax.grid(True,alpha=0.2)
    plt.savefig(plotFolder+"/scatterChargeTime"+str(suffix)+".png",transparent=False,bbox_inches='tight')
    plt.close()
    


def scatter_hist_(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')

def scatter_hist(x,y,suffix):
    fig = plt.figure(figsize=(8, 5))
    gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    # use the previously defined function
    scatter_hist_(x, y, ax, ax_histx, ax_histy)
    plt.savefig(plotFolder+"/scatterChargeTimeHist"+str(suffix)+".png",transparent=False,bbox_inches='tight')
#     plt.show()
    plt.close()

def plotScatterHist(hdfFileList,key,suffix):
    timeList,chargeList = getTimeCharge(key=key,hdf5List=hdfFileList)
    plotTimeHistogram(timeList,suffix,key)
    plotTimeChargeScatter(timeList,chargeList,suffix,key)
#     scatter_hist(timeList,chargeList,suffix)
#     print("checking length",len(timeList),len(chargeList),len(hdfFileList))


# plotScatterHist(hdf5NullList,key="OfflineIceTopSLCVEMPulses",suffix="SLCVEMTime")
# plotScatterHist(hdf5NullList,key="OfflineIceTopHLCVEMPulses",suffix="HLCVEMTime")
# plotScatterHist(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge",suffix="SLCVEMTimeCleaned")
# plotScatterHist(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge",suffix="HLCVEMTimeCleaned")
# dataT = pd.read_hdf(hdf5NullList[0],key="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration")


def getTankHitDuration(keyDuration,keyHit,hdf5List):
    hitDurationList = np.array([])
    hitList = np.array([])
    for ihdf in hdf5List:
        ihitDuration,ihit = getTankHitDuration_(keyDuration,keyHit,ihdf)
        hitDurationList = np.concatenate((hitDurationList,ihitDuration))
        hitList = np.concatenate((hitList,ihit))
    return hitDurationList,hitList

def getTankHitDuration_(keyDuration,keyHit,hdfFile):
    '''takes HLC/SLC pulses key and hdf file, returns time and charge
    '''
    dataDuration = pd.read_hdf(hdfFile,key=keyDuration)
    dataHit = pd.read_hdf(hdfFile,key=keyHit)
    return dataDuration["value"].values,dataHit["value"].values


def plotHitDuration(hit,HitDuration,prefix):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    ax.scatter(hit,HitDuration,s=10)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"tank hit", fontsize=20)
    ax.set_ylabel(r"hit duration [ns]", fontsize=20)
    # ax.set_yscale('log')
    ax.set_ylim(0.1,None)
    ax.grid(True,alpha=0.2)
    plt.savefig(plotFolder+"/scatterHitDuration"+str(prefix)+".png",transparent=False,bbox_inches='tight')
    plt.close()




# slcHitDuration,slcHit = getTankHitDuration(keyDuration="OfflineIceTopSLCTankPulsesHitTimeDuration",keyHit="OfflineIceTopSLCTankPulsesTotalHit",hdf5List=hdf5NullList)
# hlcHitDuration,hlcHit = getTankHitDuration(keyDuration="OfflineIceTopHLCTankPulsesHitTimeDuration",keyHit="OfflineIceTopHLCTankPulsesTotalHit",hdf5List=hdf5NullList)
# plotHitDuration(slcHit,slcHitDuration,prefix="SLC")
# plotHitDuration(hlcHit,hlcHitDuration,prefix="HLC")
# slcHitDuration,slcHit = getTankHitDuration(keyDuration="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration",keyHit="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalHit",hdf5List=hdf5NullList)
# hlcHitDuration,hlcHit = getTankHitDuration(keyDuration="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeHitTimeDuration",keyHit="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalHit",hdf5List=hdf5NullList)
# plotHitDuration(slcHit,slcHitDuration,prefix="SLCClean")
# plotHitDuration(hlcHit,hlcHitDuration,prefix="HLCClean")

print("slc hit w",SLCHits)
# print("slc hit w",SLCHitsW)
# plotHitDuration(SLCHitsW,SLCDuration,"SLCWeightHit")


def hitsWZenith_(hdfFile,key):
	df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
	zenithList = np.rad2deg(df_MCPrimary["zenith"].values)
	hitList = pd.read_hdf(hdfFile,key=key)["value"].values
	return hitList,zenithList

def hitsWZenith(hdfFileList,key):
	hitList = np.array([])
	zenithList = np.array([])
	for ihdf in hdfFileList:
		hits,zeniths = hitsWZenith_(ihdf,key)
		hitList = np.concatenate((hitList,hits))
		zenithList = np.concatenate((zenithList,zeniths))
	return hitList,zenithList
def plotScatterHitsZenith(hdfFileList,key,prefix):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    hitList,zenithList = hitsWZenith(hdfFileList,key)
    ax.scatter(zenithList,hitList,s=10)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
    ax.set_ylabel(r"" + str(prefix)+" hits [ns]", fontsize=20)
    # ax.set_yscale('log')
    ax.set_title(key,fontsize=16)
    ax.set_ylim(0.1,None)
    ax.grid(True,alpha=0.2)
    plt.savefig(plotFolder+"/scatterHitzenith"+str(prefix).replace(" ","")+".png",transparent=False,bbox_inches='tight')
    plt.close()

# plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopHLCVEMPulsesTotalHit",prefix="HLC VEM")
# plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopSLCVEMPulsesTotalHit",prefix="SLC VEM")

# plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalHit",prefix="Clean HLC VEM")
# plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalHit",prefix="Clean SLC VEM")

def eventsZenith(zenLim):
	dataT = pd.read_hdf(hdfFile,key="MCPrimary")
	zeniths = dataT["zenith"].values


def compareTwoWeights(hdfFileList):
	frameWeight = []
	calcWeight = weightCalc(hdfFileList)
	for ihdf in hdfFileList:
		iweightList = getWeight_(ihdf)
		frameWeight += list(iweightList)
	print("length of two weights",len(frameWeight),len(calcWeight))
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(ncols=1,nrows=1)
	ax = fig.add_subplot(gs[0])
	# hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
	hitBins = np.linspace(min(min(calcWeight),min(calcWeight)),max(max(calcWeight),max(calcWeight)),100)
	ax.hist(calcWeight,bins=hitBins,histtype="step",lw=2.5,label=r"Weight(calc)",alpha=1)
	ax.hist(frameWeight,bins=hitBins,histtype="step",lw=2.5,label="Weight(frame)",alpha=1)
	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
	ax.set_xlabel(r"weight", fontsize=20)
	ax.set_ylabel(r"count", fontsize=20)
	ax.set_yscale("log")
	# ax.set_ylim(None,10**5)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.set_title(key,fontsize=16)
	ax.legend(fontsize=10,ncol=3,loc="lower center")
	plt.savefig(plotFolder+"/weightCompare.pdf",transparent=False,bbox_inches='tight')
	plt.close()
# compareTwoWeights(hdf5NullList)

def stationRate(hdfFileList,key,prefix,hitBins,xlabel):
	stationHits = getValue(hdfFileList,key)
	weights = getValue(hdfFileList,"H4aWeight")
	triggerSTA3 = getValue(hdfFileList,"ITSMTTriggered")
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(ncols=1,nrows=1)
	ax = fig.add_subplot(gs[0])
	# hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
	ax.hist(stationHits,bins=hitBins,histtype="step",weights=[w*t for w,t in zip(weights,triggerSTA3)],lw=2.5,label=r"stations",alpha=1)
	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=12)
	ax.tick_params(axis='both',which='major',length=5)
	ax.tick_params(axis='both',which='minor',length=3)
	ax.set_xlabel(xlabel, fontsize=12)
	ax.set_ylabel(r"rate [Hz]", fontsize=12)
	ax.set_yscale("log")
	ticks = hitBins
	ax.set_xticks(ticks[::10])
	ax.set_xticks(ticks[::2],True)
	ax.set_xlim(0,None)
	# ax.set_xscale('log')
	# ax.grid(True,alpha=0.2)
	# ax.set_title(key,fontsize=16)
	# ax.legend(fontsize=10,ncol=3,loc="lower center")
	plt.savefig(plotFolder+"/Rate"+str(prefix)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# stationRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_hitStations",prefix="HLCTankStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (HLC Tank)")
# stationRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_hitStations",prefix="SLCTankStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (SLC Tank)")
# stationRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_hitTanks",prefix="HLCTankTank",hitBins=np.linspace(0,160,161),xlabel=r"$N_{tanks}$ per event (HLC Tank)")
# stationRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_hitTanks",prefix="SLCTankTank",hitBins=np.linspace(0,80,81),xlabel=r"$N_{tanks}$ per event (SLC Tank)")

# stationRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_hitStations",prefix="HLCVEMStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (HLC VEM)")
# stationRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_hitStations",prefix="SLCVEMStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (SLC VEM)")
# stationRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_hitTanks",prefix="HLCVEMTank",hitBins=np.linspace(0,160,161),xlabel=r"$N_{tanks}$ per event (HLC VEM)")
# stationRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_hitTanks",prefix="SLCVEMTank",hitBins=np.linspace(0,80,81),xlabel=r"$N_{tanks}$ per event (SLC VEM)")




def chargeRate(hdfFileList,key,prefix,hitBins,xlabel):
	Qtot = np.log10(getValue(hdfFileList,key))
	weights = getValue(hdfFileList,"H4aWeight")
	triggerSTA3 = getValue(hdfFileList,"ITSMTTriggered")
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(ncols=1,nrows=1)
	ax = fig.add_subplot(gs[0])
	# hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
	ax.hist(Qtot,bins=hitBins,histtype="step",weights=[w*t for w,t in zip(weights,triggerSTA3)],lw=2.5,label=r"stations",alpha=1)
	# ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=12)
	ax.tick_params(axis='both',which='major',length=5)
	ax.tick_params(axis='both',which='minor',length=3)
	ax.set_xlabel(xlabel, fontsize=12)
	ax.set_ylabel(r"rate [Hz]", fontsize=12)
	ax.set_yscale("log")
	ticks = hitBins
	ax.set_xticks(ticks[::10])
	ax.set_xticks(ticks[::2],True)
	# ax.set_xlim(0,None)
	# ax.set_xscale('log')
	# ax.grid(True,alpha=0.2)
	# ax.set_title(key,fontsize=16)
	# ax.legend(fontsize=10,ncol=3,loc="lower center")
	plt.savefig(plotFolder+"/chargeRate"+str(prefix)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# chargeRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalCharge",prefix="HLCTankCharge",hitBins=np.linspace(0,5,81),xlabel=r"log($Q_{tot}/VEM$) per event (HLC tank)")
# chargeRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalCharge",prefix="SLCTankCharge",hitBins=np.linspace(-2,4,81),xlabel=r"log($Q_{tot}/VEM$) per event (SLC tank)")
# chargeRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalCharge",prefix="HLCVEMCharge",hitBins=np.linspace(0,5,81),xlabel=r"log($Q_{tot}/VEM$) per event (HLC VEM)")
# chargeRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalCharge",prefix="SLCVEMCharge",hitBins=np.linspace(-2,4,81),xlabel=r"log($Q_{tot}/VEM$) per event (SLC VEM)")





