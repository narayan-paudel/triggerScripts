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
from frameTools import getGain

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
# inputList = sorted(glob.glob("/home/enpaudel/dataExp/run2023/forbushL2Window/*_Run00139197*.i3.gz"))[:1]
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

class Event(object):
  """docstring for Event"""
  def __init__(self, runID,evtID,t_start,t_end,NCh):
    super(Event, self).__init__()
    self.runID = runID
    self.evtID = evtID
    self.t_start = t_start
    self.t_end = t_end
    self.NCh = NCh

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
    self.AddParameter("GCD","GCD File","/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz")
  def Configure(self):
    self.GCD = self.GetParameter("GCD")
    f = dataio.I3File(self.GCD)
    dframe = f.pop_frame()
    while f.more() and not "I3DetectorStatus" in dframe:
      dframe = f.pop_frame()
    self.domStatus = dframe["I3DetectorStatus"].dom_status
    om = dframe["I3Geometry"].omgeo
    print("keys",[ikey for ikey in om.keys()])
    hg_doms = [ikey for ikey in om.keys() if (ikey.om in [61, 62, 63, 64] and ikey.string not in [39,74] and getGain(self.domStatus,ikey) == dataclasses.I3DOMStatus.High) or 
    (ikey.om in [62, 63, 64] and ikey.string in [39,74] and getGain(self.domStatus,ikey) == dataclasses.I3DOMStatus.High)]
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


  def nHitsHG(self,frame,pulseseries):
    NCh = 0
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
    for om,pulses in psm:
      if getGain(self.domStatus,om) == dataclasses.I3DOMStatus.High:
        NCh += len(pulses)
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
    if frame.Has("IceTopDSTPulses"):
      pulse = frame["IceTopDSTPulses"]
    else:
      print("No dst pulse in frame")
    NCh = self.nHitsHG(frame,"IceTopDSTPulses")
    thisEvent = Event(runID,evtID,start_time,end_time,NCh)
    self.evtList.append(thisEvent)

  def Finish(self):
    rateList = defaultdict(int)
    for ievt in self.evtList:
      for i,ibinstart in enumerate(self.timeBinsMJD[:-1]):
        if ievt.t_start >= self.timeBinsMJD[i] and ievt.t_start < self.timeBinsMJD[i+1]: 
          rateList[ibinstart] += ievt.NCh

    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    rateNames = ["HG DOM hits"]
    ax.errorbar(rateList.keys(),[jelt/float(self.timeBinSize) for jelt in rateList.values()],yerr=[(1.0/np.sqrt(jelt/float(self.timeBinSize))) for jelt in rateList.values()],fmt="o",ls="-",lw = 2.5,ms=8,c=colorsCustom[1],label=r"{}".format(rateNames[0]),alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.tick_params(axis='x', labelrotation=90)
    ax.set_ylabel("rate [Hz]", fontsize=20)
    ax.set_xlabel("time", fontsize=20)
    ax.set_xticks(self.timeBinsMJD,self.timeBinsHour)
    # ax.ticklabel_format(useOffset=False)
    # ax.set_yscale("log")
    # ax.set_xlim(-600,600)
    # ax.set_ylim(23.8,26) #zoom yhigh
    # ax.set_ylim(24,27) #zoom yhigh
    # ax.set_ylim(3,8) #zoom ylow
    plt.grid(visible=True, which='major', axis='both',alpha=0.6)
    ax.legend(fontsize=16,ncol=2)
    plt.savefig(plotFolder+"/domRateForbushBins.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/domRateForbushBins.png",transparent=False,bbox_inches='tight')



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
            GCD=GCD,
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

