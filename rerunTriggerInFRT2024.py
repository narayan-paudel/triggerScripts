#!/usr/bin/env python3


import os
import glob
import subprocess

from icecube.icetray import I3Tray
from icecube import icetray, dataclasses, dataio
# import icecube.icetray

from icecube import trigger_sim

from icecube import payload_parsing

import pandas as pd
import numpy as np

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default="/data/exp/IceCube/2020/filtered/level2/0508/Run00134064/Level2_IC86.2020_data_Run00134064_Subrun00000000_00000000_IT.i3.zst", help='Input files of subruns.')
args = parser.parse_args()


run2023 = "/data/sim/IceTop/2023/generated/untriggered/run2023/"

selectedFRTRunFiles = run2023+"IceTopTrig/PFFilt_PhysicsFiltering_Run00138858_Subrun00000000_00000000.i3.gz"

print("input",args.input[0].split("/")[-1])
fileName = args.input[0].split("/")[-1]


import re
run_id = [int(s) for s in re.findall(r'\d+',fileName)][0]

GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HGDomsets.i3.gz"


# outputDir = "/data/sim/IceTop/2023/generated/untriggered/run2023/FRT2024"
# outputDir = "/data/sim/IceTop/2023/generated/untriggered/run2023/FRT2023_Dec"
# outputDir = "/data/sim/IceTop/2023/generated/untriggered/run2023/FRT202X"


def checkFilter(frame,filterMask,filterKey):
    if frame.Has(filterMask):
        filterMaskObj = frame[filterMask]
        return checkFilter_(frame,filterMaskObj,filterKey)
    else:
      return False

def checkFilter_(frame,filterMaskObj,filterKey):
  if filterMaskObj[filterKey].condition_passed and filterMaskObj[filterKey].prescale_passed:
    return True
  else:
      return False

def filterTest(frame,filterKey,filterMaskKey):
  if frame.Stop == icetray.I3Frame.DAQ:
    return checkFilter(frame,filterMask=filterMaskKey,filterKey=filterKey)
  elif frame.Stop == icetray.I3Frame.Physics:
    return checkFilter(frame,filterMask=filterMaskKey,filterKey=filterKey)

class TriggerRate(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
    self.AddParameter("RunID","run id",135252)
  def Configure(self):
    self.HG7TriggerEvents = 0
    self.HLCTriggerEvents = 0
    self.nFrames = 0
    self.RunID = self.GetParameter("RunID")


  def DAQ(self,frame):
    if frame.Has("I3EventHeader"):
      self.RunID = frame["I3EventHeader"].run_id
      print("IDs",frame["I3EventHeader"].run_id,frame["I3EventHeader"].event_id)
      if frame.Has("I3TriggerHierarchy_FRT"):
        triggerHierarchy = frame['I3TriggerHierarchy_FRT']
        # print(triggerHierarchy)
        SMTTriggers = [t for t in triggerHierarchy if (t.key.config_id == 102 and t.fired)]
        HG7Triggers = [t for t in triggerHierarchy if (t.key.config_id == 30043 and t.fired)]
        if len(SMTTriggers) > 0:
          self.HLCTriggerEvents += len(SMTTriggers)
        if len(HG7Triggers) > 0:
          self.HG7TriggerEvents += len(HG7Triggers)
        self.nFrames += 1
      frame["HG7Count"] = dataclasses.I3Double(len(HG7Triggers))
      frame["HLCCount"] = dataclasses.I3Double(len(SMTTriggers))
    self.PushFrame(frame)

  def Finish(self):
    print(self.HLCTriggerEvents, self.nFrames,self.RunID)
    print("HLC6 rates {}Hz HG7 rates {}Hz {}frames".format(self.HLCTriggerEvents/self.nFrames*100,self.HG7TriggerEvents,self.nFrames))
    # with open('/data/user/enpaudel/triggerStudy/FRTRate2024/triggerRateFRT2024{}.txt'.format(self.RunID), 'a+') as f:
    with open('/data/user/enpaudel/triggerStudy/FRTRate2023/triggerRateFRT2023_{}.txt'.format(self.RunID), 'a+') as f:
    # with open('/data/user/enpaudel/triggerStudy/FRTRate202X/triggerRateFRT202X{}.txt'.format(self.RunID), 'a+') as f:
      f.write('{} {} {} {}'.format(self.RunID,self.HG7TriggerEvents,self.HLCTriggerEvents,self.nFrames))
      f.write("\n")



name = ""
QIFY = True


prekeeps = ["triggerHierarchy","OnlineFilterMask","IceTopRawData","I3Triggers","I3Triggers_FRT","I3TriggerHierarchy_FRT","I3TriggerHierarchy","I3SuperDST",
"I3EventHeader","I3DetectorStatus","I3Calibration","DSTTriggers","DOMSets"]




tray = I3Tray()
tray.AddModule("I3Reader","reader",
             filenameList=[GCD]+args.input,
            )
tray.AddModule(filterTest,"testFilter",
        # filterKey='FixedRateFilter_13',
        filterKey='FixedRateTriggerFilter',
        filterMaskKey="OnlineFilterMask",
        streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
        )

# added specially for 7HG smoke test
# decode the IceTop saved waveforms (DAQ data format in I3DAQDataIceTop
if QIFY == True:
    tray.AddModule("QConverter", "qify", WritePFrame=False)
def it_decode_needed(frame):
    if frame.Has("IceTopRawData"):
        return False
    else:
        return True
tray.AddSegment(payload_parsing.I3DOMLaunchExtractor,
                'Make_IT_launches',
                BufferID = 'I3DAQData',
                TriggerID = '',
                HeaderID = '',
                InIceID = 'UnusedInIce',
                If = it_decode_needed
                )



# tray.Add("SimpleMajorityTriggerTank","HLC6",TriggerConfigID = 102,inputHits="HLC",triggerCountName="SMT102",triggerSource=dataclasses.ICE_TOP)
# tray.Add("SimpleMajorityTriggerTank","HG7",TriggerConfigID = 30042,inputHits="both_hg",triggerCountName="IceTop7HG",triggerSource=dataclasses.ICE_TOP)
tray.AddModule("I3TriggerSimModule",
 IceTopLaunches = "IceTopRawData",
 OutputName = "I3Triggers_FRT",
 InIcePulses = "I3RecoPulseSeriesMapExtensions")

tray.AddModule("I3GlobalTriggerSim",name + "_global_trig",
  I3TriggerName="I3Triggers_FRT",
  GlobalTriggerName="I3TriggerHierarchy_FRT",
  RunID=run_id,
  FilterMode = True)
# tray.AddModule("Dump")
tray.AddModule("Keep","keep_before_write", keys = prekeeps)

tray.AddModule(TriggerRate,RunID = run_id)

# tray.AddModule("I3Writer","i3writer",
#             filename=str(outputDir)+"/"+fileName,
#             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             # DropOrphanStreams=[icetray.I3Frame.DAQ],
#             # DropOrphanStreams=[icetray.I3Frame.DAQ],
#             )
tray.Execute()
tray.Finish()
