#!/usr/bin/env python3

import os, sys
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


import healpy as hp

import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 30})

# mpl.rcParams.update({'font.size': 22})

from weighting import GetWeight, ParticleType, PDGCode
from inclinedTriggerTools import *

# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanSeedSame/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/"
basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTest/"

hdf5NullListP = sorted(glob.glob(basePath+"p*Clean*.hdf5"))
hdf5NullListHe = sorted(glob.glob(basePath+"He*Clean*.hdf5"))
hdf5NullListO = sorted(glob.glob(basePath+"O*Clean*.hdf5"))
hdf5NullListFe = sorted(glob.glob(basePath+"Fe*Clean*.hdf5"))
hdf5NullList = np.concatenate((hdf5NullListP,hdf5NullListHe,hdf5NullListO,hdf5NullListFe))
# hdf5NullList = ["/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/pDAT000001GenDetFiltProcUniqueCleanVEMEvts.hdf5"] 

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

evtList = extractEvents(hdf5NullList)
evtList_1PeV = selectEnergyEvents(evtList,10**15,0.1)
evtList_10PeV = selectEnergyEvents(evtList,10**16,0.1)

NSIDE = 32
NPIX = hp.nside2npix(NSIDE)
m = np.arange(NPIX)


def efficiencyCalc(eventList,triggerType,containment,zenith,ztol):
	"""zenith in degrees"""
	if zenith <= 65:
		if containment == True:
			# evtList = containedEvents(evtList,640)
			evtList = containedEvents(eventList,410)
		evtBin = [ievt for ievt in evtList if np.deg2rad(zenith-ztol) <= ievt.zenith < np.deg2rad(zenith+ztol)]
		totalEvts = len(evtBin)
		if triggerType == "IceTopSMT":
			triggerList = [ievt for ievt in evtBin if abs(ievt.ITSMTTriggered-1)<0.01]
		elif triggerType == "sta1":
			triggerList = [ievt for ievt in evtBin if abs(ievt.STA1Trigger-1)<0.01]
		elif triggerType == "tank3":
			triggerList = [ievt for ievt in evtBin if abs(ievt.tank3Trig-1)<0.01]
		elif triggerType == "tank4":
			triggerList = [ievt for ievt in evtBin if abs(ievt.tank4Trig-1)<0.01]
		elif triggerType == "tank5":
			triggerList = [ievt for ievt in evtBin if abs(ievt.tank5Trig-1)<0.01]
		elif triggerType == "tank6":
			triggerList = [ievt for ievt in evtBin if abs(ievt.tank6Trig-1)<0.01]
		elif triggerType == "tank7":
			triggerList = [ievt for ievt in evtBin if abs(ievt.tank7Trig-1)<0.01]
		elif triggerType == "tank8":
			triggerList = [ievt for ievt in evtBin if abs(ievt.tank8Trig-1)<0.01]
		if totalEvts == 0:
			trigEff = 0.0
		else:
			trigEff = triggerEfficiency(len(triggerList),totalEvts)
	else:
		trigEff = np.NAN
	return trigEff

def efficiencyPixels(ipix):
	"""
	ipix=pixel number
	output:zenith angel correcponding to the pixel
	"""
	theta, phi = np.degrees(hp.pix2ang(nside=NSIDE,ipix=ipix))
	return theta, phi


# print("theta",min(theta),max(theta)) #0-180
# print("phi",min(phi),max(phi)) #0-380

longitude = np.radians(np.linspace(-180, 180, 30))
latitude = np.radians(np.linspace(-90, 90, 30))

def trigEfficiencyMap(eventList,triggerType,label):
	NSIDE = 32
	NPIX = hp.nside2npix(NSIDE)
	m = np.arange(NPIX)
	theta, phi = efficiencyPixels(np.arange(NPIX))
	theta_IC = np.asarray([180-t for t in theta])
	efficiencyList = np.asarray([efficiencyCalc(eventList,triggerType=triggerType,containment=True,zenith=izenith,ztol=1) for izenith in theta_IC])
	fig = plt.figure(figsize = (10,8))
	hp.mollview(
		efficiencyList,
		# m*1.5,
		coord=["C","G"],
		unit=r"trigger efficiency",
		notext=False,
		cbar=None,
		title=r"{}".format(triggerType),
		 # norm="hist",
		 min=0,
		 max=1,
		 )
	fig = plt.gcf()
	ax = plt.gca()
	image = ax.get_images()[0]
	cmap = fig.colorbar(image, ax=ax,fraction=0.022, pad=0.005)
	cmap.ax.set_ylabel('trigger efficiency',rotation=270,fontsize=20,labelpad=23)
	cmap.ax.set_yticklabels(["{:.1f}".format(i) for i in cmap.get_ticks()])
	hp.graticule(label=True)
	plt.savefig(plotFolder+"/"+str(triggerType)+str(label)+"skymap.pdf",transparent=False,bbox_inches='tight')
	plt.close()

trigEfficiencyMap(eventList=evtList_10PeV,triggerType="IceTopSMT",label="10PeV")
trigEfficiencyMap(eventList=evtList_10PeV,triggerType="tank6",label="10PeV")
trigEfficiencyMap(eventList=evtList_1PeV,triggerType="IceTopSMT",label="1PeV")
trigEfficiencyMap(eventList=evtList_1PeV,triggerType="tank6",label="1PeV")


longitude = np.radians(np.linspace(-180, 180, 30))
latitude = np.radians(np.linspace(-90, 90, 30))



def justGrids(eventList,triggerType,label):
	NSIDE = 32
	NPIX = hp.nside2npix(NSIDE)
	m = np.arange(NPIX)
	theta, phi = efficiencyPixels(np.arange(NPIX))
	theta_IC = np.asarray([180-t for t in theta])
	efficiencyList = np.asarray([efficiencyCalc(eventList,triggerType=triggerType,containment=True,zenith=izenith,ztol=1) for izenith in longitude])
	fig = plt.figure(figsize = (15,12))
	ax = fig.add_subplot(111, projection='mollweide')
	plt.xlabel('longitude')
	plt.ylabel('latitude')
	# plt.legend()
	plt.grid(True)
	plt.savefig(plotFolder+"/TestSkyMap.pdf",transparent=False,bbox_inches='tight')
	plt.close()

# justGrids(eventList=evtList_10PeV,triggerType="IceTopSMT",label="10PeV")
