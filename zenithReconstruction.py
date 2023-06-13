#!/bin/env python3

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob

# inputFiles = ["/home/enpaudel/icecube/triggerStudy/simFiles/dataSetFullEff/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz",
# "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetFullEff/FeDAT000105GenDetFiltProcUniqueCleanVEMEvts.i3.gz"]
inputFiles = sorted(glob.glob("/home/enpaudel/icecube/triggerStudy/simFiles/dataSetFullEff/*"))
# outFile = "/home/enpaudel/icecube/triggerStudy/simFiles/wellRecoIncl.i3.gz"
outFile = "/home/enpaudel/icecube/triggerStudy/simFiles/wellRecoInclGoodCore.i3.gz"
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
zenithFullEffDict = {0.0:14.9,0.1:15.0,0.2:15.1,0.3:15.2,0.4:15.3,0.5:15.4,0.6:15.6,0.7:15.9} #zenith in rad, E in lgE(eV)

def openingAngle(theta1,phi1,theta2,phi2):
  return np.arccos(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2))

def selectWellReco(frame,zen1,zen2):
  zenith_true = frame["MCPrimary"].dir.zenith
  azimuth_true = frame["MCPrimary"].dir.azimuth
  zenith_reco = frame["ShowerPlane"].dir.zenith
  azimuth_reco = frame["ShowerPlane"].dir.azimuth
  openAngle = openingAngle(zenith_true,azimuth_true,zenith_reco,azimuth_reco)
  xcore = frame["MCPrimary"].pos.x
  ycore = frame["MCPrimary"].pos.y
  xcore_reco = frame["ShowerCOG"].pos.x
  ycore_reco = frame["ShowerCOG"].pos.y
  r_diff = np.sqrt((xcore-xcore_reco)**2+(ycore-ycore_reco)**2)/I3Units.m
  if np.arcsin(np.sqrt(zen1)) <= zenith_true < np.arcsin(np.sqrt(zen2)):
    if frame["MCPrimary"].energy*I3Units.GeV/I3Units.eV >= 10**zenithFullEffDict[zen1]:
      if np.rad2deg(zenith_true) > 60 and  np.rad2deg(openAngle) < 1 and r_diff < 20:
        print("well reconstructed inclined showers",frame["I3EventHeader"].event_id,frame["I3EventHeader"].run_id,np.rad2deg(zenith_true),np.rad2deg(openAngle),np.rad2deg(zenith_reco))
        return True
      else:
        return False
    else:
      return False
  else:
    return False


class zenithReco(icetray.I3Module):
  """docstring for zenithReco"""
  def __init__(self, ctx):
    icetray.I3Module.__init__(self,ctx)
    self.AddParameter("zenithBin","zenith Bin",[0.0,0.1])
  def Configure(self):
    self.zenithBin = self.GetParameter("zenithBin")
    self.zenithFullEffDict = {0.0:14.9,0.1:15.0,0.2:15.1,0.3:15.2,0.4:15.3,0.5:15.4,0.6:15.6,0.7:15.9} #zenith in rad, E in lgE(eV)
    self.zenithDiff = []
    self.openingAngleList = []
    self.r_diffList = []


  def Physics(self,frame):
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
    if np.arcsin(np.sqrt(self.zenithBin[0])) <= zenith_true < np.arcsin(np.sqrt(self.zenithBin[1])):
      if frame["MCPrimary"].energy*I3Units.GeV/I3Units.eV >= 10**self.zenithFullEffDict[self.zenithBin[0]]:
        # print("zenith comparision",np.arcsin(np.sqrt(self.zenithBin[0]))*180.0/np.pi,zenith_true*180.0/np.pi,np.arcsin(np.sqrt(self.zenithBin[1]))*180.0/np.pi)
        self.zenithDiff.append((zenith_true - zenith_reco)*180.0/np.pi)
        self.openingAngleList.append(openAngle*180.0/np.pi)
        if np.rad2deg(zenith_true) > 60 and  np.rad2deg(openAngle) < 1:
          print("well reconstructed inclined showers",frame["I3EventHeader"].event_id,frame["I3EventHeader"].run_id,np.rad2deg(zenith_true),np.rad2deg(openAngle),np.rad2deg(zenith_reco))
    # if np.arcsin(np.sqrt(self.zenithBin[0])) <= zenith_true < np.arcsin(np.sqrt(self.zenithBin[1])) and not np.isnan(r_diff) and np.sqrt(xcore**2+ycore**2)<410:
    # if np.arcsin(np.sqrt(self.zenithBin[0])) <= zenith_true < np.arcsin(np.sqrt(self.zenithBin[1])) and np.sqrt(xcore**2+ycore**2)<410:
    #   if frame["MCPrimary"].energy*I3Units.GeV/I3Units.eV >= 10**self.zenithFullEffDict[self.zenithBin[0]]:
    #     self.r_diffList.append(r_diff)

  def Finish(self):
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    print("zenith bin",self.zenithBin,len(self.openingAngleList))
    print("removing nan")
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    print("zenith bin",self.zenithBin,len(self.openingAngleList))
    bins = np.linspace(0,100,101)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    self.ax.hist(self.openingAngleList,bins=bins,histtype="step",label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$, E >= 10$^{{{2:.1f}}}$ eV".format(np.rad2deg(np.arcsin(np.sqrt(self.zenithBin[0]))),
      np.rad2deg(np.arcsin(np.sqrt(self.zenithBin[1]))),self.zenithFullEffDict[self.zenithBin[0]]),lw=2.5)
    self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\psi$ [$^{\circ}$]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_xlim(0,100)
    # self.ax.set_ylim(0,4300)
    # self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/zenithDiffHG7Only{:.0f}.png".format(np.rad2deg(np.arcsin(np.sqrt(self.zenithBin[0])))),transparent=False,bbox_inches='tight')
    plt.close()
    #####################################################
    #####################################################
    # self.fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=1,ncols=1)
    # self.ax = self.fig.add_subplot(gs[0])
    # # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    # # bins = np.linspace(0,90,91)
    # # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    # p68 = np.percentile(self.r_diffList,68)
    # print(p68)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$, E >= 10$^{{{2:.1f}}}$ eV".format(np.rad2deg(np.arcsin(np.sqrt(self.zenithBin[0]))),
    #   np.rad2deg(np.arcsin(np.sqrt(self.zenithBin[1]))),self.zenithFullEffDict[self.zenithBin[0]]),lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    # self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # self.ax.set_xlabel(r"abs(r$_{true}$-r$_{reco}$) [m]", fontsize=22)
    # self.ax.set_ylabel("count", fontsize=22)
    # self.ax.legend()
    # plt.savefig(plotFolder+"/coreDiffHG7Only{:.0f}.png".format(np.rad2deg(np.arcsin(np.sqrt(self.zenithBin[0])))),transparent=False,bbox_inches='tight')
    ###########################################################################
    
# for n,izen in enumerate(sin2ZenBins[:-1]):
#   tray = I3Tray()
#   tray.AddModule("I3Reader","reader",
#     filenameList=inputFiles)
#   tray.Add(zenithReco,"zenReco",
#     zenithBin=[sin2ZenBins[n],sin2ZenBins[n+1]])

#   # tray.AddModule("I3Writer","i3writer",
#   #             filename=outFile,
#   #             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#   #             )
#   tray.Execute()
#   tray.Finish()

tray = I3Tray()
tray.AddModule("I3Reader","reader",
  filenameList=inputFiles)
# tray.Add(zenithReco,"zenReco",
#   zenithBin=[0.7,0.822])
tray.Add(selectWellReco,"wellReco",
  zen1=0.7,
  zen2=0.822,
  streams=[icetray.I3Frame.Physics])

tray.AddModule("I3Writer","i3writer",
            filename=outFile,
            streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            DropOrphanStreams=[icetray.I3Frame.DAQ]
            )
tray.Execute()
tray.Finish()