#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)


import scipy.interpolate



import glob

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--vertEvents", dest="vertEvents",
                    default=False, action="store_true", required=False,
                    help="plot vertical events")
args = parser.parse_args()

# usage python hitsInclinedShower.py --vertEvents


# vertEvents = False
vertEvents = args.vertEvents
if not vertEvents:
  fileDir = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/"
  fileName = "combinedDeltaTIncl"
  plotSuffix = "Incl"
else:
  fileDir = "/home/enpaudel/dataExp/dataSetClean_VerticalLE/"
  fileName = "combinedDeltaTVert"
  plotSuffix = "Vert"

inFile = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz"
# fileList = sorted(glob.glob(fileDir+"*.i3.*"))[:2]
fileDir = "/home/enpaudel/dataExp/dataSetClean/"
fileList = sorted(glob.glob(fileDir+"*.i3.*"))[:3]

inclinationCut = 60 #degree
energyCut = 10**16 #eV
fileDir = "/home/enpaudel/dataExp/dataSetClean/"

exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

outputDir = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/"
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

def openingAngle(theta1,phi1,theta2,phi2):
  return np.arccos(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2))

def test7HG(frame):
  return frame["tank7_3000"]>0

def AddTotalCharge(frame):
  '''calculates total SLC charge in SLC tank pulses'''
  # print("Anything")
  if (frame.Has('OfflineIceTopSLCTankPulses')):
    # print("frame has OfflineIceTopSLCTankPulses")
    slc = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'OfflineIceTopSLCTankPulses')
    ITTotalChargeSLC = sum([sum([c.charge for c in slc[i]]) for i in slc.keys()])
    frame["ITTotalChargeSLC"] = dataclasses.I3Double(ITTotalChargeSLC)

def AddTotalCharge(frame,keys):
  '''calculates total SLC or HLC charges in tank pulses
  keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
  '''
  # print("Anything")
  for key in keys:    
    if (frame.Has(str(key))):
      # print("frame has", str(key))
      lc = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,str(key))
      ITTotalCharge = sum([sum([c.charge for c in lc[i]]) for i in lc.keys()])
      frame[str(key)+"TotalCharge"] = dataclasses.I3Double(ITTotalCharge)

def AddTotalTankHit(frame,pulseseriesList):
  '''calculates total SLC or HLC charges in tank pulses
  keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
  '''
  # print("Adding tank Hit")
  for pulseseries in pulseseriesList:
    NCh = 0
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap or psm.__class__ == dataclasses.I3RecoPulseSeriesMapUnion:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
    for om,pulses in psm:
      for pulse in pulses:
        NCh +=1
        break
    channels = [omkey for omkey,ps in psm if len(ps)>0]
    # if NCh == 0:
      # print("psm",psm)
    frame[str(pulseseries)+"TotalHit"] = dataclasses.I3Double(NCh)
    # print("No. of Hits in ",pulseseries,NCh,len(channels),frame[str(pulseseries)+"TotalTankHit"])

class QtotCheck(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.nHitHLCList = []
    self.nHitSLCList = []
    self.nHitList = []
    self.QtotHLCList = []
    self.QtotSLCList = []
    self.QtotList = []
    self.energyList = []
    self.zenithList = []
    self.sin2tList = []


  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    qtotHLC = frame["OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalCharge"]
    qtotSLC = frame["OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalCharge"]
    qtot = frame["IceTopTankPulsesTotalCharge"]
    nHitHLC = frame["OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalHit"]
    nHitSLC = frame["OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalHit"]
    nHit = frame["IceTopTankPulsesTotalHit"]
    # print("cores",(xcore,ycore),(xcore_reco,ycore_reco),(xcore-xcore_reco,ycore-ycore_reco),r_diff)
    zenith_true = np.rad2deg(frame["MCPrimary"].dir.zenith)
    azimuth_true = np.rad2deg(frame["MCPrimary"].dir.azimuth)
    energy_true = np.log10(frame["MCPrimary"].energy*10**9)
    # print("zenith True",zenith_reco,zenith_true,np.arcsin(np.sqrt(self.zenithBin[0])),np.arcsin(np.sqrt(self.zenithBin[1])))
    # if np.arcsin(np.sqrt(self.zenithBin[0])) <= zenith_true < np.arcsin(np.sqrt(self.zenithBin[1])) and not np.isnan(openAngle):
    self.nHitHLCList.append(nHitHLC.value)
    self.nHitSLCList.append(nHitSLC.value)
    self.nHitList.append(nHit.value)
    self.QtotHLCList.append(qtotHLC.value)
    self.QtotSLCList.append(qtotSLC.value)
    self.QtotList.append(np.log10(qtot.value))
    self.energyList.append(energy_true)
    self.zenithList.append(zenith_true)
    self.sin2tList.append(np.sin(frame["MCPrimary"].dir.zenith)**2)

  def Finish(self):
    # self.fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=1,ncols=1)
    # self.ax = self.fig.add_subplot(gs[0])
    # x = np.asarray(self.zenithList)
    # y = np.asarray(self.energyList)
    # z = np.asarray(self.QtotList)
    # # print("x",x)
    # sc = self.ax.scatter(x,y,s=z,c=z)
    # # self.ax.colorbar()
    # self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # # self.ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # self.ax.set_xlabel(r"$\theta$$^{\circ}$", fontsize=22)
    # self.ax.set_ylabel("log(E/eV)", fontsize=22)
    # cbar = self.fig.colorbar(sc, ax=self.ax)
    # cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # # cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
    # cbar.ax.set_ylabel(r"Q$_{tot}$", fontsize=22)
    # # self.ax.set_yscale("log")
    # self.ax.legend()
    # plt.savefig(plotFolder+"/qtot"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    # plt.savefig(plotFolder+"/qtot"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')

    # self.fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=1,ncols=1)
    # self.ax = self.fig.add_subplot(gs[0])
    # x = np.asarray(self.zenithList)
    # y = np.asarray(self.energyList)
    # z = np.asarray(self.QtotList)
    # # print("x",x)
    # sc = self.ax.scatter(x,y,s=z,c=z)
    # # self.ax.colorbar()
    # self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # # self.ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # self.ax.set_xlabel(r"$\theta$$^{\circ}$", fontsize=22)
    # self.ax.set_ylabel("log(E/eV)", fontsize=22)
    # cbar = self.fig.colorbar(sc, ax=self.ax)
    # cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # # cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
    # cbar.ax.set_ylabel(r"Q$_{tot}$", fontsize=22)
    # # self.ax.set_yscale("log")
    # self.ax.legend()
    # plt.savefig(plotFolder+"/qtotZen"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    # plt.savefig(plotFolder+"/qtotZen"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')

    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    x = np.asarray(self.zenithList)
    y = np.asarray(self.energyList)
    z = np.asarray(self.QtotList)
    # print("x",x)
    # sc = self.ax.scatter(x,y,s=z,c=z)
    # self.ax.colorbar()
    # self.ax.plot(self.QtotList,self.zenithList,".",ls='-',lw = 2.5,c=next(colorsIter),label=r"",alpha=1)
    # self.ax.plot(self.QtotList,self.zenithList,".",c=next(colorsIter),label=r"",alpha=1)
    self.ax.plot(self.QtotList,self.sin2tList,".",c=next(colorsIter),label=r"",alpha=1)
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # self.ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # self.ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    self.ax.set_ylabel(r"sin$^2$$\theta$", fontsize=22)
    # self.ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    self.ax.set_xlabel("log10(Q$_{tot}/[VEM])$", fontsize=22)
    # self.ax.set_yscale("log")
    self.ax.legend()
    plt.savefig(plotFolder+"/qtotZen"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/qtotZen"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')

    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    x = np.asarray(self.zenithList)
    y = np.asarray(self.energyList)
    z = np.asarray(self.QtotList)
    # print("x",x)
    # sc = self.ax.scatter(x,y,s=z,c=z)
    # self.ax.colorbar()
    # self.ax.plot(self.QtotList,self.zenithList,".",ls='-',lw = 2.5,c=next(colorsIter),label=r"",alpha=1)
    self.ax.plot(self.QtotList,self.energyList,".",c=next(colorsIter),label=r"",alpha=1)
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # self.ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # self.ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    self.ax.set_xlabel("log10(Q$_{tot}/[VEM])$", fontsize=22)
    self.ax.set_ylabel("log10(E/eV)", fontsize=22)
    # self.ax.set_yscale("log")
    self.ax.legend()
    plt.savefig(plotFolder+"/qtotE"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/qtotE"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')

    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    x = np.asarray(self.zenithList)
    y = np.asarray(self.energyList)
    z = np.asarray(self.nHitList)
    # print("x",x)
    # sc = self.ax.scatter(x,y,s=z,c=z)
    # self.ax.colorbar()
    # self.ax.plot(self.QtotList,self.zenithList,".",ls='-',lw = 2.5,c=next(colorsIter),label=r"",alpha=1)
    # self.ax.plot(self.QtotList,self.zenithList,".",c=next(colorsIter),label=r"",alpha=1)
    self.ax.plot(self.nHitList,self.sin2tList,".",c=next(colorsIter),label=r"",alpha=1)
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # self.ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # self.ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    self.ax.set_ylabel(r"sin$^2$$\theta$", fontsize=22)
    # self.ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    self.ax.set_xlabel("nHit", fontsize=22)
    # self.ax.set_yscale("log")
    self.ax.legend()
    plt.savefig(plotFolder+"/nHitZen"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/nHitZen"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')

    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    x = np.asarray(self.zenithList)
    y = np.asarray(self.energyList)
    z = np.asarray(self.QtotList)
    # print("x",x)
    # sc = self.ax.scatter(x,y,s=z,c=z)
    # self.ax.colorbar()
    # self.ax.plot(self.QtotList,self.zenithList,".",ls='-',lw = 2.5,c=next(colorsIter),label=r"",alpha=1)
    self.ax.plot(self.nHitList,self.energyList,".",c=next(colorsIter),label=r"",alpha=1)
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # self.ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # self.ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    self.ax.set_xlabel("nHit", fontsize=22)
    self.ax.set_ylabel("log10(E/eV)", fontsize=22)
    # self.ax.set_yscale("log")
    self.ax.legend()
    plt.savefig(plotFolder+"/nHitE"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/nHitE"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')




tray = I3Tray()
tray.AddModule("I3Reader","reader",
             # filenameList = args.input,
             filenameList = fileList,
             # filename = inFile,
            )

tray.AddModule(test7HG,"7HG",
              streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
              )

def Unify(frame, Keys, Output):
  """
  Simple utility to merge RecoPulseSerieses into a single Union.
  """
  extants = [k for k in Keys if k in frame]
  union = dataclasses.I3RecoPulseSeriesMapUnion(frame, extants)
  frame[Output] = union
tray.Add(Unify,"UnionHLCSLC",
  Keys=["OfflineIceTopHLCTankPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge"],
  Output='IceTopTankPulses',
  streams = [icetray.I3Frame.DAQ],
  )
tray.AddModule(AddTotalCharge,"totCharge",
               Keys=["IceTopTankPulses"],
                )
tray.AddModule(AddTotalTankHit,"totHit",
                pulseseriesList=["IceTopTankPulses"])
tray.AddModule(QtotCheck,"zTot",
            # streams = [icetray.I3Frame.DAQ],
            )

# tray.AddModule(nHitCheck,"nhits",
#             # streams = [icetray.I3Frame.DAQ],
#             )

# tray.AddModule("I3Writer","i3writer",
#             filename=str(outputDir)+fileName,
#             streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )

tray.Execute()
tray.Finish()