#!/usr/bin/env python3

import os
import glob
import subprocess

import h5py
import tables
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

import simweights

import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

from weightingPy import GetWeight, ParticleType, PDGCode
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

print("using files",hdf5NullList[:1])


# colorsList = ['#9467bd', '#e377c2','#1f77b4','#2ca02c','#bcbd22','#ff7f0e','#8c564b','#7f7f7f','#17becf','#d62728',
#       '#4477AA', '#332288', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
#       '#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']
colorsList = ['#1f77b4','#ff7f0e','#2ca02c','#8c564b','#9467bd', '#e377c2','#bcbd22','#7f7f7f','#17becf','#d62728',
      '#4477AA', '#332288','#2ca02c', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
      '#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']


colorsCustom = qualitative_colors(12)
colorsIter = iter(colorsCustom)
colorsCustom2 = colorsCustom + colorsCustom

nfilesP = nCorFiles(hdf5NullListP)
nfilesHe = nCorFiles(hdf5NullListHe)
nfilesO = nCorFiles(hdf5NullListO)
nfilesFe = nCorFiles(hdf5NullListFe)
print("no of files",nfilesP,nfilesHe,nfilesO,nfilesFe)

triggerListSelectDict = {"HG7_3000":"IT7HG","HLC6_5000":"ITSMT"}

trigWindow = 10**(-6) # in ns

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
# sin2ZenBins = [0.0,0.822]

H4aWeightList=[]
HLC6_5000List=[]

##############using simweights################

# # for hdfFile in hdf5NullList:
# #   # print(getValue_(hdfFile,key="H4aWeight"))
# #   H4aWeightList += list(getValue_(hdfFile,key="H4aWeight"))
# #   HLC6_5000List += list(getValue_(hdfFile,key="HLC6_5000"))
# # print("before weighting",sum(HLC6_5000List))
# # HLC6_5000WeightedList = []
# # for iweight,ihlc6 in zip(H4aWeightList,HLC6_5000List):
# #   HLC6_5000WeightedList.append(iweight*ihlc6)
# # print("after weighting",sum(HLC6_5000WeightedList))



# # #for use of Sim weights
# # for hdfFile in hdf5NullList:
# #   simFile = pd.HDFStore(hdfFile,"r")
# #   flux_model = simweights.GaisserH4a()
# #   weight_obj = simweights.CorsikaWeighter(simFile,100)
# #   weights = weight_obj.get_weights(flux_model)
# #   print(f"Rate = {weights.sum():5.2f} Hz")
# #   convWeights = list(getValue_(hdfFile,key="H4aWeight"))
# #   print(f"convRate = {convWeights.sum():5.2f} Hz")
# weighter = None
# hdf5NullList = ["/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/pDAT000001GenDetFiltProcUniqueCleanVEMEvtsTestWeighting.hdf5"]
# f = h5py.File(hdf5NullList[0], 'r')
# print("keys",[key for key in f.keys()])
# # for filename in hdf5NullList:
#   # filename = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/pDAT000001GenDetFiltProcUniqueCleanVEMEvtsTest.hdf5"
#   # f = h5py.File(filename, 'r')
#   # print("keys",[key for key in f.keys()])
# file_obj = tables.open_file(hdf5NullList[0], "r")
# if weighter is None:
#     weighter = simweights.IceTopWeighter(file_obj)
# else:
#     weighter += simweights.IceTopWeighter(file_obj)

# flux = simweights.GaisserH4a_IT()
# weights = weighter.get_weights(flux)

# print(weighter.tostring(flux))
# print("simweight",len(weights),sum(weights),weights)

# H4aWeightList += list(getValue_(hdf5NullList[0],key="H4aWeight"))
# print("weighting module weights",len(H4aWeightList),sum(H4aWeightList),H4aWeightList)



####################################



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
#   x,y = getCore(hdfFileList)
#   zenList,ptypeList,energyList = getZenithTypeEnergy(hdfFileList)
#   plotCoreScatter_(x,y,"all_shower","all energy")
#   energyBins = 10**(np.linspace(5.0,8.0,31))
#   print("energyBins",energyBins)
#   for n,nEnergy in enumerate(energyBins[:-1]):
#     xInBin = []
#     yInBin = []
#     for ix,iy,ienergy in zip(x,y,energyList):
#       if ienergy >= energyBins[n] and ienergy < energyBins[n+1]:
#         xInBin.append(ix) 
#         yInBin.append(iy)
#     print("nBins",n,energyBins[n],energyBins[n+1])
#     plotCoreScatter_(xInBin,yInBin,r"{:.1f}".format(np.log10(nEnergy)),r"lg(E[GeV]):{0:.1f}-{1:.1f}".format(np.log10(energyBins[n]),np.log10(energyBins[n+1])))

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
# plotScatterCore(evtList,"sta3")
# plotScatterCore(evtList,"sta1")

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
    # evtList = containedEvents(evtList,640)
    evtList = containedEvents(evtList,410)
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
    ncolor = colorsCustom2[nbin]
    for ebin, ebinStart in enumerate(energyBins[:-1]):
      lowEdge_E = energyBins[ebin]
      highEdge_E = energyBins[ebin+1]
      evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
      # totalEvts = len(evtEBin)
      weights = [ievt.H4aWeight for ievt in evtEBin]
      totalEvts = len(evtEBin)
      # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
      if triggerType == "HLC6_5000":
        triggerList = [ievt.HLC6_5000 for ievt in evtEBin]
      elif triggerType == "HG7_3000":
        triggerList = [ievt.HG7_3000 for ievt in evtEBin]
      # elif triggerType == "slc3":
      #   triggerList = [ievt.slc3Trig for ievt in evtEBin]
      # elif triggerType == "slc4":
      #   triggerList = [ievt.slc4Trig for ievt in evtEBin]
      # elif triggerType == "slc5":
      #   triggerList = [ievt.slc5Trig for ievt in evtEBin]
      trigEff = triggerEfficiency(sum(triggerList),totalEvts)
      efficiencyList.append(trigEff)
      energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))   
    ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  ax.set_ylabel(r"trigger efficiency", fontsize=22)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  # ax.set_xscale('log')
  ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  ax.set_ylim(0,1.01)
  ax.set_xlim(14,17)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  ax.grid(True,alpha=0.5)
  # ax.legend(fontsize=12)
  l1=ax.legend(loc="upper left",fontsize=12)
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='gray', marker='',markersize=5, label=r"0.98")
  l2 = ax.legend(handles=[point_dash],loc="upper right",fontsize=12,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  ax.add_artist(l1)
  ax.add_artist(l2)
  plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"Efficiency.pdf",transparent=False,bbox_inches='tight')
  plt.close()

for itrigger in triggerListSelectDict.keys():
  # plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType=itrigger,containment=False)
  plotTrigEfficiency(evtList,energyBins,weighting=False,triggerType=itrigger,containment=True)


def plotEnergyFlux(eventList,triggerType,yscale,suffix,energyScale,containment):
  """
  plots energy flux
  """
  eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-0)>0.01]
  if containment == True:
    # eventList = containedEvents(eventList,410)
    eventList = containedEvents(eventList,410)
    # eventList = containedEvents(eventList,800)
  # eventList = [ievt for ievt in eventList if 10**6.9 <= ievt.energy < 10**8.0]    
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(ncols=1,nrows=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # hitBins = np.linspace(14.0,17.0,31)
  hitBins = np.linspace(14.0,16.8,29)
  totalRate = 0
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
    energyList = []
    neventList = []
    # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
    ncolor = next(colorIter)
    ax,histSum = plotSteps(evtZenBin,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),energyScale,ncolor=colorsCustom2[nbin])
    totalRate += histSum
  # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
  # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
  ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
  ax.set_yscale(yscale) 
  ax.set_ylim(10**-5,10**0.0)
  # ax.set_xlim(14.0,17.0)
  ax.set_xlim(14.0,16.8)
  # ax.set_xscale('log')
  ax.grid(True,alpha=0.7)
  ax.text(0.82,0.82,s=r"{0} rate:{1:.1f} Hz".format(triggerListSelectDict[triggerType],totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
  # ax.set_title(key,fontsize=16)
  ax.legend(fontsize=10,ncol=3,loc="lower center")
  plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+".png",transparent=False,bbox_inches='tight')
  plt.close()

def plotSteps(triggeredEvts,ax,legendLabel,energyScale,ncolor):
  hitBins = np.linspace(14.0,17.0,31)
  energy = [ievt.energy for ievt in triggeredEvts]
  # weights = [ievt.H4aWeight for ievt in triggeredEvts]
  weights = [ievt.simWeight for ievt in triggeredEvts]
  # weights_direct = [ievt.directWeight for ievt in triggeredEvts]
  energy = np.log10(energy)+9.0
  hist,binEdge = np.histogram(energy,hitBins,weights=[w for w in weights])
  # binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
  print("sum of simulated rate",sum(hist))
  # if "total" in legendLabel:
    # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
  # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
    # ax.set_ylim(10**-4,20)
  binCenter = (binEdge[:-1]+binEdge[1:])/2.0
  if str(energyScale) == "0.0":
    H = [h for h in hist]   
  else:
    H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
  # ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
  # ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.6f} Hz".format(sum(hist)),color=ncolor,alpha=1)
  ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.3f} Hz".format(sum(hist)),color=ncolor,alpha=1)
  return ax,sum(hist)


for itrigger in triggerListSelectDict.keys():
  plotEnergyFlux(evtList,triggerType=itrigger,yscale="log",suffix="fluxLog",energyScale=0.0,containment=False)
  plotEnergyFlux(evtList,triggerType=itrigger,yscale="log",suffix="fluxLog",energyScale=0.0,containment=True)
