#!/usr/bin/env python3

from icecube.icetray import I3Tray
from icecube import icetray,dataclasses,dataio

import re
import glob
import subprocess

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input',"-i",nargs="+",type=str,
  default="/data/exp/IceCube/2023/filtered/PFFilt/1231/PFFilt_PhysicsFiltering_Run00138807_Subrun00000000_00000043.tar.bz2",
  help="input data GCD for IceTop")
parser.add_argument('--GCD',"-g",type=str,
  default="/home/enpaudel/dataExp/run2023/GCD/PFGCD_Run00138807_Subrun00000000.i3.gz",help="GCD file")
args = parser.parse_args()

import numpy as np

GCD = args.GCD
inputList = args.input

# fileName = PFFilt_PhysicsFiltering_Run00138937_Subrun00000000_00000120.tar.bz2
# sub_run = re.findall(r'\d+',args.input[0].split("/")[-1])[2]
# sub_run = re.findall(r'\d+',args.input[0].split("/")[-1])
sub_run = re.findall(r'\d+',args.input[0].split("/")[-1])[2]
# sub_run = args.input[0].split("Subrun00000000_")[-1].split(".")[0]
# print("sub_run",sub_run)
filename = args.input[0].split("/")[-1]
print("filename",filename)
# outputDir = "/data/sim/IceTop/2023/generated/untriggered/run2023/IceTopTrig/"
outputDir = "/data/sim/IceTop/2023/generated/untriggered/run2023/forbush/"

class TriggerRate(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.HG7TriggerEvents = 0
    self.HG7_HLCTriggerEvents = 0
    self.HG7_NoHLCTriggerEvents = 0
    self.HLCTriggerEvents = 0


  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    if frame.Has("I3EventHeader"):
      self.runID = frame["I3EventHeader"].run_id
      if frame.Has("DSTTriggers"):
        triggerHierarchy = frame['DSTTriggers'].unpack(frame['I3DetectorStatus'])
        SMTTriggers = [t for t in triggerHierarchy if (t.key.config_id == 102 and t.fired)]
        HG7Triggers = [t for t in triggerHierarchy if (t.key.config_id == 30043 and t.fired)]
        if len(SMTTriggers) > 0:
          self.HLCTriggerEvents += 1
        if len(HG7Triggers) > 0:
          self.HG7TriggerEvents += 1
          if len(SMTTriggers) < 1:
            self.HG7_NoHLCTriggerEvents += 1
          elif len(SMTTriggers) > 0:
            self.HG7_HLCTriggerEvents += 1

  def Finish(self):
    with open('/data/user/enpaudel/triggerStudy/rateFilesEvents/triggerEventRate2024Test{}.txt'.format(self.runID), 'a+') as f:
    # with open('/data/user/enpaudel/triggerStudy/triggerEventRateForbush{}.txt'.format(self.runID), 'a+') as f:
      f.write('{} {} {} {} {} {}'.format(self.runID,sub_run,self.HG7TriggerEvents,self.HLCTriggerEvents,self.HG7_NoHLCTriggerEvents,self.HG7_HLCTriggerEvents))
      f.write("\n")

def selectIceTopTrigger(frame):
  hasIT = False
  if frame.Has("DSTTriggers"):
    triggerHierarchy = frame['DSTTriggers'].unpack(frame['I3DetectorStatus'])
    frame["triggerHierarchy"] = triggerHierarchy
    if len(triggerHierarchy) > 0:
      icetop_triggers = [t for t in triggerHierarchy if t.key.source == dataclasses.ICE_TOP and t.fired]
      if len(icetop_triggers) > 0:
        hasIT = True
  return hasIT




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
