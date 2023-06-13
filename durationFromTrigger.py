#!/usr/bin/env python3


import os
import sys
from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube import topeventcleaning, tpx,toprec

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from customColors import qualitative_colors

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+',
 default="/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTanks/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz",
 help='Input files after running detector.py.')
args = parser.parse_args()

GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"

exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

colorsCustom = qualitative_colors(12)

#example python durationHits.py /home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTanks/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz
#example python durationFromTrigger.py /home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTanks/FeDAT005957GenDetFiltProcUniqueCleanVEMEvts.i3.gz
#example python durationFromTrigger.py /home/enpaudel/icecube/triggerStudy/background/i3_daq/hit7HG/tank7EventsPulse/Tank7Events_0000.i3.gz
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

class checkDeltaT(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.SMTTriggered = 0
    self.SMTUntriggered = 0
    self.GlobalTriggered = 0
    self.GlobalUntriggered = 0
    self.totalEvents = 0
    self.delta_t_lowEVert = []
    self.delta_t_lowEIncl = []
    self.delta_t_highEVert = []
    self.delta_t_highEIncl = []


  def checkTriggers(self,frame,triggerHierarchy,config_id):
    trigLabel = ConfigIDMap[config_id][2]+"_"+str(ConfigIDMap[config_id][1])
    if len(triggerHierarchy) == 0:
      frame[trigLabel] = dataclasses.I3Double(0)
    else:
      SMTTriggers = [t for t in triggerHierarchy if (t.key.config_id == config_id and t.fired)]
      if len(SMTTriggers) != 0:
        frame[trigLabel] = dataclasses.I3Double(1)
      else:
        frame[trigLabel] = dataclasses.I3Double(0)
    return frame

  def deltaTHitFromTrigger(self,frame,pulseseries):
    # print(frame.keys())
    if frame.Has("I3TriggerHierarchy"):
      triggerHierarchy = frame["I3TriggerHierarchy"]
    elif frame.Has("QTriggerHierarchy"):
      triggerHierarchy = frame["QTriggerHierarchy"]
    # TankTriggers = [t for t in triggerHierarchy if (t.key.config_id == 30043 and t.fired)]
    TankTriggers = [t for t in triggerHierarchy if (t.key.config_id == 124 and t.fired)]
    if len(TankTriggers) == 0:
      TankTriggers = [t for t in triggerHierarchy if (t.key.config_id == 30043 and t.fired)]
    # print("TankTriggers",len(TankTriggers),TankTriggers)
    if len(TankTriggers) == 1:
      triggerTime = TankTriggers[0].time
    else:
      print("There are multiple {} triggers".format(TankTriggers[0].key.config_id))
    # print("TankTriggers",[t.key.config_id for t in triggerHierarchy])
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap or dataclasses.I3RecoPulseSeriesMapUnion:
      # print(frame[pulseseries].sources)
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
      # psm = frame[pulseseries].apply(frame)
    times_list = []
    hit_stations = []
    hit_omkeys = []
    for om,pulses in psm:
      hit_stations.append(om[0])
      hit_omkeys.append(om)
      # break
    hit_stations = list(set(hit_stations))
    hit_tanks = []
    hit_times = []
    pulsesList = []
    for istation in hit_stations:
      tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
      (iom[0] == istation and iom[0] in [26,39,74] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
      pulsesList += [ipulse for om in tanks for ipulse in psm[om]]
    pulseList = sorted(pulsesList, key=lambda pulse: pulse.time)
    hit_times = [pulse.time for pulse in pulseList]
    shifted_hit_times = []
    if len(hit_times) > 0:
      shifted_hit_times = [itime-triggerTime for itime in hit_times]
    # print("trigger time",triggerTime,hit_times,shifted_hit_times)
    return shifted_hit_times

  def deltaTHit(self,frame,pulseseries):
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap or dataclasses.I3RecoPulseSeriesMapUnion:
      # print(frame[pulseseries].sources)
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
      # psm = frame[pulseseries].apply(frame)
    times_list = []
    hit_stations = []
    hit_omkeys = []
    for om,pulses in psm:
      hit_stations.append(om[0])
      hit_omkeys.append(om)
      # break
    hit_stations = list(set(hit_stations))
    hit_tanks = []
    hit_times = []
    pulsesList = []
    for istation in hit_stations:
      tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
      (iom[0] == istation and iom[0] in [26,39,74] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
      pulsesList += [ipulse for om in tanks for ipulse in psm[om]]
    pulseList = sorted(pulsesList, key=lambda pulse: pulse.time)
    hit_times = [pulse.time for pulse in pulseList]
    deltaT = np.array([])
    if len(hit_times) > 0:
      deltaT = np.diff(hit_times)
    return list(deltaT)

  # def DAQ(self,frame):
  #   """what if there is only one trigger in trigger hierarchy
  #   """
  #   # dT = self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
  #   if frame["SMT273"] == 1:
  #     dT = self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
  #     if np.rad2deg(frame["MCPrimary"].dir.zenith) <= 45 and 10**6.2<=frame["MCPrimary"].energy<10**6.8:
  #       self.delta_t_lowEVert += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
  #     if np.rad2deg(frame["MCPrimary"].dir.zenith) > 45 and 10**6.2<=frame["MCPrimary"].energy<10**6.8:
  #       self.delta_t_lowEIncl += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
  #     if np.rad2deg(frame["MCPrimary"].dir.zenith) <= 45 and 10**7.0<=frame["MCPrimary"].energy<10**7.5:
  #       self.delta_t_highEVert += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
  #     if np.rad2deg(frame["MCPrimary"].dir.zenith) > 45 and 10**7.0<=frame["MCPrimary"].energy<10**7.5:
  #       self.delta_t_highEIncl += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
  #   # print(frame.keys())
  #   self.totalEvents += 1
  #   self.PushFrame(frame)

  def DAQ(self,frame):
    """what if there is only one trigger in trigger hierarchy
    """
    # dT = self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
    if frame.Has("SMT273"):      
      if frame["SMT273"] == 1:
        dT = self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
        if np.rad2deg(frame["MCPrimary"].dir.zenith) <= 45 and 10**6.2<=frame["MCPrimary"].energy<10**6.8:
          self.delta_t_lowEVert += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
        if np.rad2deg(frame["MCPrimary"].dir.zenith) > 45 and 10**6.2<=frame["MCPrimary"].energy<10**6.8:
          self.delta_t_lowEIncl += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
        if np.rad2deg(frame["MCPrimary"].dir.zenith) <= 45 and 10**7.0<=frame["MCPrimary"].energy<10**7.5:
          self.delta_t_highEVert += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
        if np.rad2deg(frame["MCPrimary"].dir.zenith) > 45 and 10**7.0<=frame["MCPrimary"].energy<10**7.5:
          self.delta_t_highEIncl += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
    else:

      TankTriggers = [t for t in frame["I3TriggerHierarchy"] if (t.key.config_id == 30043 and t.fired)]
      if len(TankTriggers) == 1:
        dT = self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
        self.delta_t_lowEVert += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
        self.delta_t_lowEIncl += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
        self.delta_t_highEVert += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
        self.delta_t_highEIncl += self.deltaTHitFromTrigger(frame,"IceTopVEMPulses")
    # print(frame.keys())
    self.totalEvents += 1
    self.PushFrame(frame)


  def getLog(self,x):
    # return [np.log10(i) for i in x if i!=0]
    return x

  def Finish(self):
    print("total number of events",self.totalEvents)
    colorIter = iter(colorsCustom)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(-1,12000,12002)
    # bins = None
    # bins = np.linspace(np.log10(1),np.log10(12000),30)
    # bins = np.linspace(-1,2000,2002)
    # bins = np.linspace(-1000,12000,13002)
    # bins = np.linspace(-100,3500,180)
    # bins = np.linspace(-100,12000,500)
    bins = np.linspace(-1*10**4,1*10**4,5*10**2)
    density = False
    # bins = np.linspace(0,10000,501)
    ax.hist(self.getLog(self.delta_t_lowEVert),bins=bins,density=density,histtype="step",color=next(colorIter),label=r"15.2$\leq$lgE[eV]<15.8, $0^{\circ}<\theta\leq 45^{\circ}$",lw=1.5)
    ax.hist(self.getLog(self.delta_t_lowEIncl),bins=bins,density=density,histtype="step",color=next(colorIter),label=r"15.2$\leq$lgE[eV]<15.8, $45^{\circ}<\theta\leq 65^{\circ}$",lw=1.5)
    ax.hist(self.getLog(self.delta_t_highEVert),bins=bins,density=density,histtype="step",color=next(colorIter),label=r"16.0$\leq$lgE[eV]<16.5, $0^{\circ}<\theta\leq 45^{\circ}$",lw=1.5)
    ax.hist(self.getLog(self.delta_t_highEIncl),bins=bins,density=density,histtype="step",color=next(colorIter),label=r"16.0$\leq$lgE[eV]<16.5, $45^{\circ}<\theta\leq 65^{\circ}$",lw=1.5)
    # ax.hist(self.delta_t_lowEVert,bins=bins,density=density,histtype="step",color=next(colorIter),label=r"15.2$\leq$lgE[eV]<15.8, $0^{\circ}<\theta\leq 45^{\circ}$",lw=1.5)
    # ax.hist(self.delta_t_lowEIncl,bins=bins,density=density,histtype="step",color=next(colorIter),label=r"15.2$\leq$lgE[eV]<15.8, $45^{\circ}<\theta\leq 65^{\circ}$",lw=1.5)
    # ax.hist(self.delta_t_highEVert,bins=bins,density=density,histtype="step",color=next(colorIter),label=r"16.0$\leq$lgE[eV]<16.5, $0^{\circ}<\theta\leq 45^{\circ}$",lw=1.5)
    # ax.hist(self.delta_t_highEIncl,bins=bins,density=density,histtype="step",color=next(colorIter),label=r"16.0$\leq$lgE[eV]<16.5, $45^{\circ}<\theta\leq 65^{\circ}$",lw=1.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # ax.set_xlabel("delta_t", fontsize=22)
    ax.set_xlabel("t_hit-t_trigger[ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    ax.set_yscale("log")
    # ax.set_xlim()
    ax.grid(True,alpha=0.4)
    # pulseDict = {"nSLC":"SLC tank","nHLC":"HLC tank","nHLCVEM":"HLC VEM","nSLCVEM":"SLC VEM"}
    # ax.text(0.82,0.88,s="mean {0} hits\n{1:.1f} per event".format(pulseDict[hitType],meanHit),size=15,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    # ax.legend(loc="upper left",fontsize=12)
    ax.legend(fontsize=8)
    # ax.hist(weightspy3,histtype="step")
    plt.savefig(plotFolder+"/durationFromTrigger.pdf",transparent=False,bbox_inches='tight')
    plt.close()


tray = I3Tray()
tray.AddModule("I3Reader","reader",
  filenameList=[GCD]+args.input,
            )
keep_list = ['I3DetectorStatus','I3Geometry','I3Triggers','CleanIceTopRawData','H4aWeight',
    'HLC6_5000','tank7_3000','I3EventHeader','QTriggerHierarchy','QFilterMask',"I3TriggerHierarchy"
  'IceTopRawData', 'MCPrimary','IceTopPulses',
  'OfflineIceTopHLCVEMPulses', 'OfflineIceTopHLCTankPulses','OfflineIceTopSLCTankPulses', 'OfflineIceTopSLCVEMPulses',
   "SMT183", "SMT263", "SMT273", "SMT283",
   "DrivingTime","ExcludedHLCTanks","TankKey","ExcludedSLCTanks",
   "TankKey","FADC","I3EventHeader","I3EventHeader","I3TriggerHierarchy",
   "I3Trigger","ITHLCPulseInfo","ITHLCTankPulses","ITHLCVEMPulses",
   "ITHLCWaveforms","ITSLC","ITSLCTankPulses","ITSLCVEMPulses",
   "ITWaveformRange","I3TimeWindow","ITWaveforms","IceTopHLCPEPulses",
   "IceTopRawData","IceTopSLCPEPulses","IceTopTankPulses",
   "I3RecoPulseSeriesMapUnion","IceTopVEMPulses",
   "I3RecoPulseSeriesMapUnion","InIceRawData"
   ]

tray.AddModule("Keep","keep unusual time",
             keys = keep_list,
             # If = True,
             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
             )
def Unify(frame, Keys, Output):
  """
  Simple utility to merge RecoPulseSerieses into a single Union.
  """
  extants = [k for k in Keys if k in frame]
  union = dataclasses.I3RecoPulseSeriesMapUnion(frame, extants)
  frame[Output] = union

# tray.Add(Unify,"UnionHLCSLC",
#   Keys=['OfflineIceTopHLCTankPulses','OfflineIceTopSLCTankPulses'],
#   Output='IceTopTankPulses',
#   streams=[icetray.I3Frame.DAQ],
#   )

# tray.Add(Unify,"UnionHLCSLCVEM",
#   Keys=['OfflineIceTopHLCVEMPulses','OfflineIceTopSLCVEMPulses'],
#   Output='IceTopVEMPulses',
#   streams=[icetray.I3Frame.DAQ],
#   )
tray.AddModule(checkDeltaT, "dT",
        # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
        )



# tray.AddModule("I3Writer","i3writer",
#   filename=str(outputDir)+dataSetClean+str(fileName)+"CleanVEMEvts.i3.gz",
#   streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )
tray.Execute()
tray.Finish()