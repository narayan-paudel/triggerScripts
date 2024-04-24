#!/usr/bin/env python3

from icecube.icetray import I3Tray
from icecube import icetray,dataclasses,dataio

import re
import glob
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import datetime

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('--input',"-i",nargs="+",type=str,
#   default="/data/exp/IceCube/2023/filtered/PFFilt/1231/PFFilt_PhysicsFiltering_Run00138807_Subrun00000000_00000043.tar.bz2",
#   help="input data GCD for IceTop")
# parser.add_argument('--GCD',"-g",type=str,
#   default="/home/enpaudel/dataExp/run2023/GCD/PFGCD_Run00138807_Subrun00000000.i3.gz",help="GCD file")
# args = parser.parse_args()

import numpy as np

# GCD = args.GCD
GCD = "/home/enpaudel/dataExp/run2023/GCD/PFGCD_Run00139197_Subrun00000000.i3.gz"
# inputList = args.input
# inputList = sorted(glob.glob("/home/enpaudel/dataExp/run2023/forbush/*.i3.gz"))[:1]
inputList = sorted(glob.glob("/home/enpaudel/dataExp/run2023/forbushL2/*Unpack.i3.gz"))
# inputList = [irun for  irun in inputList if "139197" in irun][:5]
print("inputList",inputList)
runDur = {139196:8*60*60+0*60+11.164081,139197:8*60*60+0*60+10.310869,139198:8*60*60+0*60+11.468356}
dates = [np.datetime64("2023-11-28T15:11:43"),np.datetime64("2024-03-24T23:11:55"),np.datetime64("2024-03-25T07:12:05.45")]

# fileName = PFFilt_PhysicsFiltering_Run00138937_Subrun00000000_00000120.tar.bz2
# sub_run = re.findall(r'\d+',args.input[0].split("/")[-1])[2]
# sub_run = re.findall(r'\d+',args.input[0].split("/")[-1])
sub_run = re.findall(r'\d+',inputList[0].split("/")[-1])[2]
# sub_run = args.input[0].split("Subrun00000000_")[-1].split(".")[0]
# print("sub_run",sub_run)
filename = inputList[0].split("/")[-1]
print("filename",filename)
# outputDir = "/data/sim/IceTop/2023/generated/untriggered/run2023/IceTopTrig/"
# outputDir = "/data/sim/IceTop/2023/generated/untriggered/run2023/forbush/"

class TriggerRate(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.HG7Triggers = {}
    self.SMTTriggers = {}
    self.HG7TriggersHLC6 = {}
    self.SMTTriggersHLC6 = {}
    self.HG7TriggersHLC4 = {}
    self.SMTTriggersHLC4 = {}
    self.HG7TriggersHLC2 = {}
    self.SMTTriggersHLC2 = {}
    self.HG7TriggersHLC0 = {}
    self.SMTTriggersHLC0 = {}
    self.dictList = [self.HG7Triggers,self.HG7TriggersHLC6,self.HG7TriggersHLC4,self.HG7TriggersHLC2,self.HG7TriggersHLC0,self.SMTTriggers,self.SMTTriggersHLC6,self.SMTTriggersHLC4,self.SMTTriggersHLC2,self.SMTTriggersHLC0]
    self.dictListName = ["HG7Triggers","HG7TriggersHLC6","HG7TriggersHLC4","HG7TriggersHLC2","HG7TriggersHLC0","SMTTriggers","SMTTriggersHLC6","SMTTriggersHLC4","SMTTriggersHLC2","SMTTriggersHLC0"]
    # self.dictList = [self.HG7Triggers,self.HG7TriggersHLC6,self.HG7TriggersHLC4,self.HG7TriggersHLC2]
    # self.dictListName = ["HG7Triggers","HG7TriggersHLC6","HG7TriggersHLC4","HG7TriggersHLC2"]

    for irun in [139196,139197,139198]:
      for idict in self.dictList:
        idict[irun] = 0
  def nHits(self,frame,pulseseries):
    NCh = 0
    psm = frame["IceTopHLCVEMPulses"]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, "IceTopHLCVEMPulses")
    for om,pulses in psm:
      for pulse in pulses:
        NCh +=1
        break
    channels = [omkey for omkey,ps in psm if len(ps)>0]
    # print(NCh, len(channels))
    return NCh



  # def Physics(self,frame):
  def DAQ(self,frame):
    runID = frame["I3EventHeader"].run_id
    if frame.Has("IceTopHLCVEMPulses"):
      VEMPulse_HLC = frame["IceTopHLCVEMPulses"]
    else:
      print("No HLC VEM pulse in frame")
    NCh_HLC = self.nHits(frame,"IceTopHLCVEMPulses")
    if frame.Has("IceTopSLCVEMPulses"):
      VEMPulse_HLC = frame["IceTopHLCVEMPulses"]
    else:
      print("No SLC VEM pulse in frame")
    NCh_SLC = self.nHits(frame,"IceTopSLCVEMPulses")

    # if int(frame["I3EventHeader"].event_id) == int(8351):
    if frame.Has("I3EventHeader"):
      runID = frame["I3EventHeader"].run_id
      print(runID)
      if frame.Has("DSTTriggers"):
        for trigger in frame['DSTTriggers'].unpack(frame['I3DetectorStatus']):
            if trigger.key.config_id == 102 and trigger.fired:
                self.SMTTriggers[runID] += 1
                if NCh_HLC < 1:
                  self.SMTTriggersHLC0[runID] += 1
                if NCh_HLC == 2:
                  self.SMTTriggersHLC2[runID] += 1
                if NCh_HLC == 4:
                  self.SMTTriggersHLC4[runID] += 1
                if NCh_HLC >= 6:
                  self.SMTTriggersHLC6[runID] += 1
            if trigger.key.config_id == 30043 and trigger.fired:
                self.HG7Triggers[runID] += 1
                if NCh_HLC < 1:
                  self.HG7TriggersHLC0[runID] += 1
                if NCh_HLC == 2:
                  self.HG7TriggersHLC2[runID] += 1
                if NCh_HLC == 4:
                  self.HG7TriggersHLC4[runID] += 1
                if NCh_HLC >= 6:
                  self.HG7TriggersHLC6[runID] += 1
  def Finish(self):
    print(self.HG7Triggers)
    runList = [139196,139197,139198]
    for runID in runList:      
      print("{}  HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7Triggers[runID]/runDur[runID],self.SMTTriggers[runID]/runDur[runID]))
      print("{} 6 HLC HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7TriggersHLC6[runID]/runDur[runID],self.SMTTriggersHLC6[runID]/runDur[runID]))
      print("{} 4 HLC HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7TriggersHLC4[runID]/runDur[runID],self.SMTTriggersHLC4[runID]/runDur[runID]))
      print("{} 2 HLC HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7TriggersHLC2[runID]/runDur[runID],self.SMTTriggersHLC2[runID]/runDur[runID]))
      print("{} 0 HLC HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7TriggersHLC0[runID]/runDur[runID],self.SMTTriggersHLC0[runID]/runDur[runID]))
      with open('/data/user/enpaudel/triggerStudy/triggerRateForbush{}.txt'.format(runID), 'a+') as f:
        f.write('{} {} {} {} {} {} {} {}'.format(runID,self.HG7Triggers[runID],self.SMTTriggers[runID],self.HG7TriggersHLC6[runID],self.SMTTriggersHLC6[runID],self.HG7TriggersHLC4[runID],self.SMTTriggersHLC4[runID],self.HG7TriggersHLC2[runID],self.SMTTriggersHLC2[runID]))
        f.write("\n")

    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for n,idict in enumerate(self.dictList):
      HG7rateList = []
      for irun in runList:
        HG7rateList.append(idict[runID]/runDur[runID])
      ax.plot(runList,HG7rateList,".",label=self.dictListName[n],alpha=1)
    # ax.plot(runList,HLCrateList,".",label=r"ITSMT",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    secax = ax.secondary_xaxis('top')
    secax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    # ax.tick_params(axis='x', labelrotation=90)
    # secax.tick_params(axis='x', labelrotation=90)
    ax.set_xticks(runList[::])
    secax.set_xticks(runList[::])
    # ax.xaxis.get_major_locator().set_params(integer=True)
    secax.set_xticklabels(dates)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("rate [Hz]", fontsize=11)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel("Run", fontsize=11)
    secax.set_xlabel("day", fontsize=11)
    ax.ticklabel_format(useOffset=False)
    # ax.set_yscale("log")
    # ax.set_xlim(-600,600)
    # ax.set_ylim(-600,600)
    ax.legend()
    plt.savefig(plotFolder+"/trigCountForbush.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/trigCountForbush.png",transparent=False,bbox_inches='tight')






# def selectIceTopTrigger(frame):
#   hasIT = False
#   if frame.Has("DSTTriggers"):
#     triggerHierarchy = frame['DSTTriggers'].unpack(frame['I3DetectorStatus'])
#     frame["triggerHierarchy"] = triggerHierarchy
#     if len(triggerHierarchy) > 0:
#       icetop_triggers = [t for t in triggerHierarchy if t.key.source == dataclasses.ICE_TOP and t.fired]
#       if len(icetop_triggers) > 0:
#         hasIT = True
#   return hasIT




tray = I3Tray()
tray.AddModule("I3Reader","reader",
              FilenameList=[GCD]+inputList,
              # filenameList=inputList[0],
              # filename=GCD,
              )
tray.AddModule(TriggerRate, "rT",
            # GCD=GCD,
            # Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
            )
# tray.AddModule(selectIceTopTrigger, "ITTrig",
#             # GCD=GCD,
#             Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
#             )

# tray.AddModule("I3Writer","i3writer",
#             filename=str(outputDir)+filename.split("/")[-1].split(".")[0] + ".i3.gz",
#             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )

tray.Execute()
tray.Finish()
