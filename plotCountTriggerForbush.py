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
import math

from astropy.time import Time, TimeDelta

from collections import defaultdict

from customColors import qualitative_colors

colorsCustom = qualitative_colors(5)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
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
# inputList = sorted(glob.glob("/home/enpaudel/dataExp/run2023/forbushL2/*Unpack.i3.gz"))[:1]
inputList = sorted(glob.glob("/home/enpaudel/dataExp/run2023/forbushL2Window/*_Run00139197*.i3.gz"))
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
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
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
      VEMPulse_SLC = frame["IceTopSLCVEMPulses"]
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

class Event(object):
  """docstring for Event"""
  def __init__(self, runID,evtID,t_start,t_end,nHLC,hasHG7,hasSMT):
    super(Event, self).__init__()
    self.runID = runID
    self.evtID = evtID
    self.t_start = t_start
    self.t_end = t_end
    self.nHLC = nHLC
    self.hasHG7 = hasHG7
    self.hasSMT = hasSMT


def getTime(t_day):
  hr = t_day*24.0
  hr_int = math.trunc(hr)
  minute = (hr -hr_int)*60.0
  minute_int = math.trunc(minute)
  seconds = (minute - minute_int)*60.0
  print(hr_int,minute_int,seconds)
  return hr_int,minute_int,seconds

def getTimeSecond(t_seconds):
  # hr = t_seconds/(60.*60.)
  # hr_int = math.trunc(hr)
  # minute = (hr -hr_int)*60.0
  # minute_int = math.trunc(minute)
  # seconds = (minute - minute_int)*60.0
  # print(hr_int,minute_int,seconds)
  minutes, seconds = divmod(t_seconds, 60)
  hr_int, minute_int = divmod(minutes, 60)
  return int(hr_int%24),int(minute_int),seconds


class TriggerRateBins(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.timeBinSize = 600 #seconds
    nBins = int(10800/self.timeBinSize)
    # self.timeBins = np.linspace(0,10800,181) #3 hours binned in 60 seconds each
    self.timeBins = np.linspace(0,10800,nBins+1) #3 hours binned in 600 seconds each
    print("timebins",self.timeBins)
    # t_ref = Time("2024-03-24 23:00:00").mjd
    t_ref = Time("2024-03-24 23:50:00").mjd #this time goes -->
    self.timeBinsMJD = [t_ref + (ielt/(86400.)) for ielt in self.timeBins] #3 hours
    timeBinsSeconds = [ielt+(23*60+50)*60 for ielt in self.timeBins] # <--
    print("timeBinsSeconds",timeBinsSeconds)
    self.timeBinsHour = [r"{:02d}:{:02d}:{:03.1f}".format(*getTimeSecond(ielt)) for ielt in timeBinsSeconds]
    # self.timeBinsHour = [TimeDelta(ielt) for ielt in timeBinsSeconds]
    print(self.timeBinsMJD)
    print(self.timeBinsHour)
    self.evtList = []

  def nHits(self,frame,pulseseries):
    NCh = 0
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
    for om,pulses in psm:
      for pulse in pulses:
        NCh +=1
        break
    channels = [omkey for omkey,ps in psm if len(ps)>0]
    # print(NCh, len(channels))
    return NCh

  # def Physics(self,frame):
  def DAQ(self,frame):
    hasHG7 = False
    hasSMT = False
    evh = frame["I3EventHeader"]
    start_time = evh.start_time.mod_julian_day+evh.start_time.mod_julian_sec/(24.*3600.)+(evh.start_time.mod_julian_nano_sec*math.pow(10,-9))/(24.*3600.)
    end_time = evh.end_time.mod_julian_day+evh.end_time.mod_julian_sec/(24.*3600.)+(evh.end_time.mod_julian_nano_sec*math.pow(10,-9))/(24.*3600.)
    runID = evh.run_id
    evtID = evh.event_id
    # print("astropyTime",Time("2024-03-24 23:00:00").mjd)    
    if frame.Has("IceTopHLCVEMPulses"):
      VEMPulse_HLC = frame["IceTopHLCVEMPulses"]
    else:
      print("No HLC VEM pulse in frame")
    NCh_HLC = self.nHits(frame,"IceTopHLCVEMPulses")
    if frame.Has("triggerHierarchy"):
      trigger = frame["triggerHierarchy"]
    elif frame.Has("DSTTriggers"):
      trigger = frame['DSTTriggers'].unpack(frame['I3DetectorStatus'])
    else:
      print("missing trigger hierarchy")

    for itrigger in trigger:
      if itrigger.key.config_id == 102 and itrigger.fired:
        hasSMT = True
      if itrigger.key.config_id == 30043 and itrigger.fired:
        hasHG7 = True
    thisEvent = Event(runID,evtID,start_time,end_time,NCh_HLC,hasHG7,hasSMT)
    self.evtList.append(thisEvent)

  def Finish(self):
    # print(self.HG7Triggers)
    # runList = [139196,139197,139198]
    # for runID in runList:      
    #   print("{}  HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7Triggers[runID]/runDur[runID],self.SMTTriggers[runID]/runDur[runID]))
    #   print("{} 6 HLC HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7TriggersHLC6[runID]/runDur[runID],self.SMTTriggersHLC6[runID]/runDur[runID]))
    #   print("{} 4 HLC HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7TriggersHLC4[runID]/runDur[runID],self.SMTTriggersHLC4[runID]/runDur[runID]))
    #   print("{} 2 HLC HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7TriggersHLC2[runID]/runDur[runID],self.SMTTriggersHLC2[runID]/runDur[runID]))
    #   print("{} 0 HLC HG7 {:.1f} Hz, SMT {:.1f} Hz".format(runID,self.HG7TriggersHLC0[runID]/runDur[runID],self.SMTTriggersHLC0[runID]/runDur[runID]))
    #   with open('/data/user/enpaudel/triggerStudy/triggerRateForbush{}.txt'.format(runID), 'a+') as f:
    #     f.write('{} {} {} {} {} {} {} {}'.format(runID,self.HG7Triggers[runID],self.SMTTriggers[runID],self.HG7TriggersHLC6[runID],self.SMTTriggersHLC6[runID],self.HG7TriggersHLC4[runID],self.SMTTriggersHLC4[runID],self.HG7TriggersHLC2[runID],self.SMTTriggersHLC2[runID]))
    #     f.write("\n")
    # print("evtList",len(self.evtList))
    # print("evtList",[ievt.t_start for ievt in self.evtList])
    # print("timebins",self.timeBins)
    rateList = defaultdict(int)
    rateListHLC0 = defaultdict(int)
    rateListHLC2 = defaultdict(int)
    rateListHLC4 = defaultdict(int)
    rateListHLC6 = defaultdict(int)
    for ievt in self.evtList:
      if ievt.hasHG7:
        for i,ibinstart in enumerate(self.timeBinsMJD[:-1]):
          if ievt.t_start >= self.timeBinsMJD[i] and ievt.t_start < self.timeBinsMJD[i+1]: 
            rateList[ibinstart] += 1
            if ievt.nHLC == 0:
              rateListHLC0[ibinstart] += 1
            elif ievt.nHLC == 2:
              rateListHLC2[ibinstart] += 1
            elif ievt.nHLC == 4:
              rateListHLC4[ibinstart] += 1
            elif ievt.nHLC >= 6:
              rateListHLC6[ibinstart] += 1
    #   ax.plot(runList,HG7rateList,".",label=self.dictListName[n],alpha=1)   
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    rateNames = ["7HG"]
    for iname,ielt in enumerate([rateList]):
      # ax.plot(ielt.keys(),[jelt/float(self.timeBinSize) for jelt in ielt.values()],ls='-', marker='.',lw=2.5,ms=8,label=r"{}".format(rateNames[iname]),alpha=1)
      ax.errorbar(ielt.keys(),[jelt/float(self.timeBinSize) for jelt in ielt.values()],yerr=[(1.0/np.sqrt(jelt/float(self.timeBinSize))) for jelt in ielt.values()],fmt="o",ls="-",lw = 2.5,ms=8,c=colorsCustom[iname],label=r"{}".format(rateNames[iname]),alpha=1)
      print("ratenames",rateNames[iname])
      print("data",[i for i in zip(ielt.keys(),[jelt/float(self.timeBinSize) for jelt in ielt.values()])])
      # ax.plot(self.timeBinsHour,[jelt/60.0 for jelt in ielt.values()],".",label=r"{}".format(rateNames[iname]),alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    # secax = ax.secondary_xaxis('top')
    # secax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.tick_params(axis='x', labelrotation=90)
    # secax.tick_params(axis='x', labelrotation=90)
    # ax.xaxis.get_major_locator().set_params(integer=True)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("rate [Hz]", fontsize=20)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel("time", fontsize=20)
    # secax.set_xlabel("day", fontsize=11)
    # ax.set_xticklabels(self.timeBinsHour)
    ax.set_xticks(self.timeBinsMJD,self.timeBinsHour)
    # ax.ticklabel_format(useOffset=False)
    # ax.set_yscale("log")
    # ax.set_xlim(-600,600)
    # ax.set_ylim(23.8,26) #zoom yhigh
    ax.set_ylim(23,28) #zoom yhigh
    # ax.set_ylim(3,6.2) #zoom ylow
    plt.grid(visible=True, which='major', axis='both',alpha=0.6)
    ax.legend(fontsize=18,ncol=1)
    plt.savefig(plotFolder+"/trigCountForbushBins7HG.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/trigCountForbushBins7HG.png",transparent=False,bbox_inches='tight')

    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    rateNames = ["7HG HLC0","7HG HLC2","7HG HLC4","7HG HLC6+"]
    for iname,ielt in enumerate([rateListHLC0,rateListHLC2,rateListHLC4,rateListHLC6]):
      # ax.plot(ielt.keys(),[jelt/float(self.timeBinSize) for jelt in ielt.values()],ls='-', marker='.',lw=2.5,ms=8,label=r"{}".format(rateNames[iname]),alpha=1)
      ax.errorbar(ielt.keys(),[jelt/float(self.timeBinSize) for jelt in ielt.values()],yerr=[(1.0/np.sqrt(jelt/float(self.timeBinSize))) for jelt in ielt.values()],fmt="o",ls="-",lw = 2.5,ms=8,c=colorsCustom[iname+1],label=r"{}".format(rateNames[iname]),alpha=1)
      print("ratenames",rateNames[iname])
      print("data",[i for i in zip(ielt.keys(),[jelt/float(self.timeBinSize) for jelt in ielt.values()])])
      # ax.plot(self.timeBinsHour,[jelt/60.0 for jelt in ielt.values()],".",label=r"{}".format(rateNames[iname]),alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    # secax = ax.secondary_xaxis('top')
    # secax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.tick_params(axis='x', labelrotation=90)
    # secax.tick_params(axis='x', labelrotation=90)
    # ax.xaxis.get_major_locator().set_params(integer=True)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("rate [Hz]", fontsize=20)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel("time", fontsize=20)
    # secax.set_xlabel("day", fontsize=11)
    # ax.set_xticklabels(self.timeBinsHour)
    ax.set_xticks(self.timeBinsMJD,self.timeBinsHour)
    # ax.ticklabel_format(useOffset=False)
    # ax.set_yscale("log")
    # ax.set_xlim(-600,600)
    # ax.set_ylim(23.8,26) #zoom yhigh
    # ax.set_ylim(24,27) #zoom yhigh
    ax.set_ylim(3,8) #zoom ylow
    plt.grid(visible=True, which='major', axis='both',alpha=0.6)
    ax.legend(fontsize=16,ncol=2)
    plt.savefig(plotFolder+"/trigCountForbushBinsNHLC.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/trigCountForbushBinsNHLC.png",transparent=False,bbox_inches='tight')




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
# tray.AddModule(TriggerRate, "rT",
#             # GCD=GCD,
#             # Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
#             )
tray.AddModule(TriggerRateBins, "rT",
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
