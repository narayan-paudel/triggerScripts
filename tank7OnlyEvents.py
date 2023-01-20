#!/usr/bin/env python3

import os
import glob
import subprocess
import re

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

from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units

from customColors import qualitative_colors

import numpy as np

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

from weighting import GetWeight, ParticleType, PDGCode
from inclinedTriggerTools import *

from icecube.weighting.fluxes import GaisserH4a_IT


#!/usr/bin/env python3


# gcd_file = "/home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz"
gcd_file = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"

for f in dataio.I3File(gcd_file,'r'):
  if f.Stop == icetray.I3Frame.Geometry:
    geom = f['I3Geometry']
    # print("geometry",geom.keys())
    stageo = geom.stationgeo
    omgeo = geom.omgeo
    tankgeo = geom.tankgeo
    tank_x = [[],[]]
    tank_y = [[],[]]
    tank_z = [[],[]]
    stations = []
    omkeys = []
    for stnkey,station in geom.stationgeo:
      stations.append(stnkey)
      # print("station key",stnkey,station)
      for tank in station:
        for omkey in tank.omkey_list:
          if omkey[1] == 61:
            tank_x[0].append(tank.position.x)
            tank_y[0].append(tank.position.y)
            tank_z[0].append(tank.position.z)
          elif omkey[1] == 63:
            tank_x[1].append(tank.position.x)
            tank_y[1].append(tank.position.y)
            tank_z[1].append(tank.position.z)

def getCircle(radius):
  theta = np.linspace(0,2*np.pi,100)
  return radius*np.cos(theta),radius*np.sin(theta)

def tankLayout(ax,x,y,stations,r):
  """additionally adds a circle of radius r m """
  ax.plot(x[0],y[0],"o",ms=3,mew=1,mfc='none',c=qualitative_colors(5)[1],label="tank A",alpha=1)
  ax.plot(x[1],y[1],"o",ms=3,mew=1,mfc='none',c=qualitative_colors(5)[2],label="tank B",alpha=1)
  xCirc,yCirc = getCircle(r)
  ax.plot(xCirc,yCirc,'-',c=qualitative_colors(5)[3],lw=3.0,label="r = {:.1f} m".format(r))
  ax.set_xlabel(r"x [m]", fontsize=20)
  ax.set_ylabel(r"y [m]", fontsize=20)
  ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
  ax.tick_params(which='both', width=1.5)
  ax.tick_params(which='major', length=7)
  ax.tick_params(which='minor', length=4)
  # ax.set_yscale('log')
  ax.grid(True,alpha=0.5)
  ax.set_aspect("equal")
  ax.set_ylim(-650,650)
  ax.set_xlim(-650,650)
  # for ix,iy,ista in zip(x[0],y[0],stations):
  #   ax.text(ix-40,iy-40,s=r"{}".format(ista),size=12)
  # ax.legend(fontsize=14,ncol=2)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  return ax



# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanSeedSame/"
basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanFRT/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTest/"

hdf5NullListP = sorted(glob.glob(basePath+"p*Clean*.hdf5"))
hdf5NullListHe = sorted(glob.glob(basePath+"He*Clean*.hdf5"))
hdf5NullListO = sorted(glob.glob(basePath+"O*Clean*.hdf5"))
hdf5NullListFe = sorted(glob.glob(basePath+"Fe*Clean*.hdf5"))
# hdf5NullList = np.concatenate((hdf5NullListP,hdf5NullListHe,hdf5NullListO,hdf5NullListFe))
# hdf5NullList = np.concatenate((hdf5NullListP,hdf5NullListHe,hdf5NullListO,hdf5NullListFe))
hdf5NullList = hdf5NullListFe
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

# colorsList = ['#9467bd', '#e377c2','#1f77b4','#2ca02c','#bcbd22','#ff7f0e','#8c564b','#7f7f7f','#17becf','#d62728',
#       '#4477AA', '#332288', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
#       '#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']
colorsList = ['#1f77b4','#ff7f0e','#2ca02c','#8c564b','#9467bd', '#e377c2','#bcbd22','#7f7f7f','#17becf','#d62728',
      '#4477AA', '#332288','#2ca02c', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
      '#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray',
      '#1f77b4','#ff7f0e','#2ca02c','#8c564b','#9467bd', '#e377c2','#bcbd22','#7f7f7f','#17becf','#d62728',
      '#4477AA', '#332288','#2ca02c', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
      '#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)

# multiple_color = mpl.cm.tab20(range(20))
# multiple_color = np.concatenate((np.asarray([]),mpl.cm.winter(range(4)),mpl.cm.cool(range(4)),mpl.cm.Wistia(range(4)),mpl.cm.copper(range(4))))
multiple_color = ["olive","darkslategray","teal","darkturquoise","paleturquoise","midnightblue","royalblue","cornflowerblue","lightsteelblue","brown","indianred","lightcoral","rosybrown","darkgoldenrod","goldenrod","gold","khaki","purple","mediumorchid","violet","plum"]

# triggerList = ["HLC6_5000","tank6_5000","tank6_4000","tank6_3000","tank6_2000",
#   "tank7_5000","tank7_4000","tank7_3000","tank7_2000","tank8_5000","tank8_4000","tank8_3000","tank8_2000"]
triggerList2 = ["HLC6_5000","tank7_3000"]
# triggerList7 = ["HLC6_5000","tank7_5000","tank7_4000","tank7_3000","tank7_2000"]
# triggerListSelect = ["HLC6_5000","tank6_3000","tank7_3000","tank8_3000"]

print("proton")
nfilesP = nCorFiles(hdf5NullListP)
print("helium")
nfilesHe = nCorFiles(hdf5NullListHe)
print("Oxygen")
nfilesO = nCorFiles(hdf5NullListO)
print("Fe")
nfilesFe = nCorFiles(hdf5NullListFe)
print("no of files p He O Fe",nfilesP,nfilesHe,nfilesO,nfilesFe)
# print("no of files  He",nfilesHe)

trigWindow = 10**(-6) # in ns

# sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
sin2ZenBins = [0.7,0.822]
# sin2ZenBins = [0.0,0.822]

evtList = extractEvents(hdf5NullList)
evtList_contained = containedEvents(evtList,410)
total_events = sum([1 for ievt in evtList])
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

vertEvents = [ievt.eventID for ievt in evtList_contained if (abs(getattr(ievt,"HLC6_5000")-1)<0.01 and 0.0 <= ievt.zenith < 0.2 and abs(np.log10(getattr(ievt,"energy"))-7.0)<=0.05)]
inclEvents = [ievt.eventID for ievt in evtList_contained if (abs(getattr(ievt,"tank7_3000")-1)<0.01 and 0.95 <= ievt.zenith < 1.0 and abs(np.log10(getattr(ievt,"energy"))-7.0)<=0.05)]
inclEventsNoHLC6 = [ievt.eventID for ievt in evtList_contained if (abs(getattr(ievt,"tank7_3000")-1)<0.01 and abs(getattr(ievt,"HLC6_5000")-0)<0.01 and 0.95 <= ievt.zenith < 1.0 and abs(np.log10(getattr(ievt,"energy"))-7.0)<=0.05)]

print("vertical events",vertEvents)
print("incl events",inclEvents)
print("incl events No HLC",inclEventsNoHLC6)



# energyBins = 10**np.linspace(5, 8.0, 31)
# energyBinsShort = 10**np.linspace(6.1, 6.5, 7)
energyBinsShort = [10**7.1, 10**7.5]
# energyBins = 10**np.linspace(5, 8.0, 7)
energyBinslgE = np.linspace(5.0,8.9,4000)
energyBinCenter = [5.1,6.1,7.1,8.1]

def plotHist(x,x1,x2,xlabel,plotlabel):
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  # ax.hist(weightsi3,histtype="step")
  ax.hist(x,histtype="step",label="total",lw=2.5)
  ax.hist(x1,histtype="step",label="HLC6",lw=2.5)
  ax.hist(x2,histtype="step",label="tank7",lw=2.5)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(xlabel, fontsize=22)
  ax.set_ylabel("count", fontsize=22)
  # ax.legend(loc="upper left",fontsize=12)
  ax.legend(fontsize=12)
  # ax.hist(weightspy3,histtype="step")
  plt.savefig(plotFolder+"/hits"+str(plotlabel)+".pdf",transparent=False,bbox_inches='tight')
  plt.close()

def plotScatter(x,x1,x2,xlabel):
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  # ax.hist(weightsi3,histtype="step")
  ax = tankLayout(ax,tank_x,tank_y,stations,410)
  # ax.plot([ix[0] for ix in x],[ix[1] for ix in x],ls="",marker="o",ms=3,mew=1,mfc='none',label="total",alpha=0.6)  
  ax.plot([ix[0] for ix in x2],[ix[1] for ix in x2],ls="",marker="D",ms=3,mew=1,mfc='none',color="orange",label="tank7",alpha=0.6)
  ax.plot([ix[0] for ix in x1],[ix[1] for ix in x1],ls="",marker="*",ms=3,mew=1,mfc='none',color="plum",label="hlc6",alpha=0.6)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  # ax.set_xlabel(x, fontsize=22)
  # ax.set_ylabel(y, fontsize=22)
  # ax.legend(loc="upper right",fontsize=12)
  ax.legend(ncol=3,fontsize=8)
  # ax.hist(weightspy3,histtype="step")
  plt.savefig(plotFolder+"/Scatter"+str(xlabel)+".pdf",transparent=False,bbox_inches='tight')
  plt.close()

def plotScatterBoth(x,x1,x2,xlabel):
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  # ax.hist(weightsi3,histtype="step")
  ax = tankLayout(ax,tank_x,tank_y,stations,410)
  # ax.plot([ix[0] for ix in x],[ix[1] for ix in x],ls="",marker="o",ms=3,mew=1,mfc='none',label="total",alpha=0.6)  
  ax.plot([ix[0] for ix in x if (ix in x1 and ix in x2)],[ix[1] for ix in x if (ix in x1 and ix in x2)],ls="",marker="p",ms=2.5,mew=0.5,mfc='none',color="orange",label="both",alpha=1)
  ax.plot([ix[0] for ix in x if (ix not in x1 and ix in x2)],[ix[1] for ix in x if (ix not in x1 and ix in x2)],ls="",marker="*",ms=2.5,mew=0.5,mfc='none',color="plum",label="tank7 only",alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  # ax.set_xlabel(x, fontsize=22)
  # ax.set_ylabel(y, fontsize=22)
  # ax.legend(loc="upper right",fontsize=12)
  ax.legend(ncol=3,fontsize=8)
  # ax.hist(weightspy3,histtype="step")
  plt.savefig(plotFolder+"/Scatter"+str(xlabel)+"OnlyBoth.pdf",transparent=False,bbox_inches='tight')
  plt.close()




def plotTank7Events(evtList,energyBins,triggerType1,triggerType2,containment):
  '''
  plots trigger efficiency in different zenith bins
  '''
  if containment == True:
    # evtList = containedEvents(evtList,640)
    evtList = containedEvents(evtList,410)
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # colorIter = iter(colorsList)
  lowEdge = np.arcsin(np.sqrt(sin2ZenBins[0]))
  highEdge = np.arcsin(np.sqrt(sin2ZenBins[1]))
  evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
  energyList = []
  efficiencyList = []
  # ncolor = colorsCustom2[nbin]
  lowEdge_E = energyBins[0]
  highEdge_E = energyBins[1]
  evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
  # totalEvts = len(evtEBin)
  weights = [ievt.H4aWeight for ievt in evtEBin]
  totalEvts = len(evtEBin)
  print("total events",totalEvts)
  # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
  triggerList1 = [ievt for ievt in evtEBin if abs(getattr(ievt,triggerType1)-1)<0.01]
  triggerList2 = [ievt for ievt in evtEBin if abs(getattr(ievt,triggerType2)-1)<0.01]
  trigEff1 = triggerEfficiency(len(triggerList1),totalEvts)
  trigEff2 = triggerEfficiency(len(triggerList2),totalEvts)
  # print("triggerEfficiency",triggerType1,len(triggerList1),trigEff1,triggerType2,len(triggerList2),trigEff2)
  # print("eventID runID",[(ievt.runID,ievt.eventID,ievt.CRType) for ievt in evtEBin if (abs(getattr(ievt,triggerType2)-1)<0.01 and abs(getattr(ievt,triggerType1)-0)<0.01 and abs(np.log10(getattr(ievt,"energy"))-7.45)<=0.05)])
  # print("eventID runID both triggers",[(ievt.runID,ievt.eventID,ievt.CRType) for ievt in evtEBin if (abs(getattr(ievt,triggerType2)-1)<0.01 and abs(getattr(ievt,triggerType1)-1)<0.01 and abs(np.log10(getattr(ievt,"energy"))-7.45)<=0.05)])
  # print("eventID runID both triggers tank 7",[(ievt.runID,ievt.eventID,ievt.CRType) for ievt in evtEBin if (abs(getattr(ievt,triggerType2)-1)<0.01 and abs(np.log10(getattr(ievt,"energy"))-7.45)<=0.05)])
  print("eventID runID both triggers HLC6",[(ievt.runID,ievt.eventID,ievt.CRType) for ievt in evtEBin if (abs(getattr(ievt,triggerType1)-1)<0.01 and abs(np.log10(getattr(ievt,"energy"))-7.45)<=0.05)])
  print("eventID runID both triggers HLC6  and tank7",[(ievt.runID,ievt.eventID,ievt.CRType) for ievt in evtEBin if (abs(getattr(ievt,triggerType2)-1)<0.01 and abs(getattr(ievt,triggerType1)-1)<0.01 and abs(np.log10(getattr(ievt,"energy"))-7.45)<=0.05)])
  print("eventID runID both triggers HLC6  or tank7",[(ievt.runID,ievt.eventID,ievt.CRType) for ievt in evtEBin if ((abs(getattr(ievt,triggerType2)-1)<0.01 or abs(getattr(ievt,triggerType1)-1)<0.01) and abs(np.log10(getattr(ievt,"energy"))-7.45)<=0.05)])
  # print("crtype",[ievt.CRType for ievt in evtEBin if (abs(getattr(ievt,triggerType2)-1)<0.01 and abs(getattr(ievt,triggerType1)-0)<0.01)])
  # print("energy",[ievt.energy for ievt in evtEBin if (abs(getattr(ievt,triggerType2)-1)<0.01 and abs(getattr(ievt,triggerType1)-0)<0.01)])
  # print("zenith",[ievt.zenith for ievt in evtEBin if (abs(getattr(ievt,triggerType2)-1)<0.01 and abs(getattr(ievt,triggerType1)-0)<0.01)])
  # print("core",[(ievt.coreX,ievt.coreY) for ievt in evtEBin if (abs(getattr(ievt,triggerType2)-1)<0.01 and abs(getattr(ievt,triggerType1)-0)<0.01)])
  zenithTotal = [ievt.zenith*180.0/np.pi for ievt in evtEBin]
  zenithHLC6 = [ievt.zenith*180.0/np.pi for ievt in triggerList1]
  zenithTank7 = [ievt.zenith*180.0/np.pi for ievt in triggerList2]
  plotHist(zenithTotal,zenithHLC6,zenithTank7,r"zenith [$^{\circ}$]","zenith")
  energyTotal = [np.log10(ievt.energy*10**9) for ievt in evtEBin]
  energyHLC6 = [np.log10(ievt.energy*10**9) for ievt in triggerList1]
  energyTank7 = [np.log10(ievt.energy*10**9) for ievt in triggerList2]
  plotHist(energyTotal,energyHLC6,energyTank7,"log10 (energy [eV])","energy")
  coreTotal = [(ievt.coreX,ievt.coreY) for ievt in evtEBin]
  coreHLC6 = [(ievt.coreX,ievt.coreY) for ievt in triggerList1]
  coreTank7 = [(ievt.coreX,ievt.coreY) for ievt in triggerList2]
  plotScatter(coreTotal,coreHLC6,coreTank7,"core")
  plotScatterBoth(coreTotal,coreHLC6,coreTank7,"core")
  # ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
  # ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  # ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  # ax.set_ylabel(r"trigger efficiency", fontsize=22)
  # # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  # ax.text(0.78,0.1,s=r"trig:{0}".format(triggerType),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # # ax.set_xscale('log')
  # ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  # ax.set_ylim(0,1.01)
  # ax.set_xlim(14,17)
  # # ax.yaxis.set_minor_locator(MultipleLocator(100))
  # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  # ax.grid(True,alpha=0.5)
  # l1=ax.legend(loc="upper left",fontsize=12)
  # point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"0.98")
  # # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  # # ax.add_artist(l1)
  # # ax.add_artist(l2)
  # plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"EfficiencyTest.pdf",transparent=False,bbox_inches='tight')
  plt.close()


plotTank7Events(evtList,energyBinsShort,triggerType1="HLC6_5000",triggerType2="tank7_3000",containment=True)