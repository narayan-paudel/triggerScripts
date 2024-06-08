#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})



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
fileList = sorted(glob.glob(fileDir+"*.i3.*"))

inclinationCut = 60 #degree
energyCut = 10**16 #eV
fileDir = "/data/user/enpaudel/triggerStudy/simFiles/dataSetClean1_6/"

exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

# outputDir = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/"
outputDir = "/home/enpaudel/dataExp/dataSetClean_InclinedHE_7HG_reco/"
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

def openingAngle(theta1,phi1,theta2,phi2):
  return np.arccos(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2))

def test7HG(frame):
  return frame["HG7_3"]>0

def excludeITSMT(frame):
  return frame["HLC6"]<1



class zenithCheck(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.openingAngleList = []
    self.r_diffList = []


  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    xcore = frame["MCPrimary"].pos.x
    ycore = frame["MCPrimary"].pos.y
    xcore_reco = frame["ShowerCOG"].pos.x
    ycore_reco = frame["ShowerCOG"].pos.y
    r_diff = np.sqrt((xcore-xcore_reco)**2+(ycore-ycore_reco)**2)/I3Units.m
    # print("cores",(xcore,ycore),(xcore_reco,ycore_reco),(xcore-xcore_reco,ycore-ycore_reco),r_diff)
    zenith_true = frame["MCPrimary"].dir.zenith
    azimuth_true = frame["MCPrimary"].dir.azimuth
    zenith_reco = frame["ShowerPlane"].dir.zenith
    azimuth_reco = frame["ShowerPlane"].dir.azimuth
    openAngle = openingAngle(zenith_true,azimuth_true,zenith_reco,azimuth_reco)
    # print("zenith True",zenith_reco,zenith_true,np.arcsin(np.sqrt(self.zenithBin[0])),np.arcsin(np.sqrt(self.zenithBin[1])))
    # if np.arcsin(np.sqrt(self.zenithBin[0])) <= zenith_true < np.arcsin(np.sqrt(self.zenithBin[1])) and not np.isnan(openAngle):
    self.openingAngleList.append(openAngle*180.0/np.pi)
    self.r_diffList.append(r_diff)

  def Finish(self):
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    bins = np.linspace(-1,100,102)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    self.ax.hist(self.openingAngleList,bins=bins,histtype="step",label=r"",lw=2.5)
    self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}$^{{\circ}}$".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\psi$ [$^{\circ}$]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_xlim(0,100)
    self.ax.set_ylim(0.9,5*10**3)
    self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/openAngleHG7Only"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/openAngleHG7Only"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')
    plt.close()
    #####################################################
    #####################################################
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    p68 = np.percentile(self.r_diffList,68)
    print(p68)
    bins = np.linspace(-1,2000,2002)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"abs(r$_{true}$-r$_{reco}$) [m]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_yscale("log")
    self.ax.legend()
    plt.savefig(plotFolder+"/coreDiffHG7Only"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreDiffHG7Only"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')






tray = I3Tray()
tray.AddModule("I3Reader","reader",
             # filenameList = args.input,
             filenameList = fileList,
             # filename = inFile,
            )

tray.AddModule(test7HG,"7HG",
              streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
              )

tray.AddModule(excludeITSMT,"notITSMT",
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
tray.AddModule(zenithCheck,"zhits",
            # streams = [icetray.I3Frame.DAQ],
            )

tray.AddModule("I3Writer","i3writer",
            filename=str(outputDir)+fileName+".i3.gz",
            streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            )

tray.Execute()
tray.Finish()