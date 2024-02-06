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
else:
  fileDir = "/home/enpaudel/dataExp/dataSetClean_VerticalLE/"
  fileName = "combinedDeltaTVert"

inFile = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz"
fileList = sorted(glob.glob(fileDir+"*.i3.*"))

inclinationCut = 60 #degree
energyCut = 10**16 #eV
fileDir = "/home/enpaudel/dataExp/dataSetClean/"

exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

outputDir = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/"
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

def test7HG(frame):
  print(frame["tank7_3000"])
  print("test",frame["tank7_3000"]>0)
  return frame["tank7_3000"]>0



class timeCheck(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
    self.AddParameter("pulseseriesList","ps",["IceTopTankPulses"])
  def Configure(self):
    self.pulseseriesList = self.GetParameter("pulseseriesList")
    self.hitTimes_diff = np.array([])
    self.hitTimes_duration = []


  def DAQ(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    hitTimesEvent = np.array([])
    for pulseseries in self.pulseseriesList:
      psm = frame[pulseseries]
      if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap or psm.__class__ == dataclasses.I3RecoPulseSeriesMapUnion:
        psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
      hit_stations = []
      hit_omkeys = []
      for om,pulses in psm:
        for pulse in pulses:
          hit_stations.append(om[0])
          hit_omkeys.append(om)
          break
      hit_stations = list(set(hit_stations))
      hit_tanks = []
      for istation in hit_stations:
        tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
        (iom[0] == istation and iom[0] in [26,39,74] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
        pulses = [psm[om] for om in tanks]
        hit_tanks += tanks
        # ipulses = [ipulse for ipulse in pulse for pulse in pulses]
        hit_times = np.sort([pulse[0].time for pulse in pulses])
        # hit_times = [ipulse.time for pulse in pulses for ipulse in pulse]
        hitTimesEvent = np.concatenate((hitTimesEvent,hit_times))
    hitTimesEvent = np.sort(hitTimesEvent)
    if len(hitTimesEvent) >= 2:
      deltaT = np.diff(hitTimesEvent)
      self.hitTimes_diff = np.concatenate((self.hitTimes_diff,deltaT))
      self.hitTimes_duration.append(hitTimesEvent[-1]-hitTimesEvent[0])

  def Finish(self):
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # bins = np.linspace(0,10000,10001)
    bins = np.linspace(-1,2000,2002)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    self.ax.hist(self.hitTimes_diff,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.hitTimes_diff,histtype="step",label=r"",lw=2.5)
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # self.ax.set_xlabel(r"$\psi$ [$^{\circ}$]", fontsize=22)
    self.ax.set_xlabel(r"time [ns]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    # self.ax.set_xlim(0,100)
    # self.ax.set_ylim(10**0,3*10**4)
    self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/"+fileName+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/"+fileName+".pdf",transparent=False,bbox_inches='tight')
    plt.close()
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # bins = np.linspace(0,10000,10001)
    bins = np.linspace(-1,15000,15002)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    self.ax.hist(self.hitTimes_duration,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.hitTimes_duration,histtype="step",label=r"",lw=2.5)
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # self.ax.set_xlabel(r"$\psi$ [$^{\circ}$]", fontsize=22)
    self.ax.set_xlabel(r"time [ns]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    # self.ax.set_xlim(0,100)
    # self.ax.set_ylim(10**0,3*10**4)
    self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/"+fileName+"Dur.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/"+fileName+"Dur.pdf",transparent=False,bbox_inches='tight')
    plt.close()






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
tray.AddModule(timeCheck,"t_hits",
            pulseseriesList=["IceTopTankPulses"],
            # streams = [icetray.I3Frame.DAQ],
            )

# tray.AddModule("I3Writer","i3writer",
#             filename=str(outputDir)+fileName,
#             streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )

tray.Execute()
tray.Finish()