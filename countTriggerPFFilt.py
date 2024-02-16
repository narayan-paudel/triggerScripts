#!/usr/bin/env python3

from icecube.icetray import I3Tray
from icecube import icetray,dataclasses,dataio

import re
import glob
import subprocess

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input',"-i",nargs="+",type=str,
  default="/home/enpaudel/dataExp/run2023/IT7HGTrig/PFFilt_PhysicsFiltering_Run00138616_Subrun00000000_00000197_7HG.i3.gz",
  help="input data GCD for IceTop")
parser.add_argument('--GCD',"-g",type=str,
  default="/home/enpaudel/dataExp/run2023/GCD/PFGCD_Run00138616_Subrun00000000.i3.gz",help="GCD file")
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

class TriggerRate(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.HG7Triggers = 0
    self.HLCTriggers = 0

  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    if frame.Has("I3EventHeader"):
      self.runID = frame["I3EventHeader"].run_id
      if frame.Has("DSTTriggers"):
        for trigger in frame['DSTTriggers'].unpack(frame['I3DetectorStatus']):
            # if trigger.key.config_id == 30043 and trigger.fired:
            if trigger.key.config_id == 102 and trigger.fired:
                self.HLCTriggers += 1
            if trigger.key.config_id == 30043 and trigger.fired:
                self.HG7Triggers += 1
  def Finish(self):
    with open('/data/user/enpaudel/triggerStudy/triggerRate{}.txt'.format(self.runID), 'a+') as f:
      f.write('{} {} {} {}'.format(self.runID,sub_run,self.HG7Triggers,self.HLCTriggers))
      f.write("\n")



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
tray.Execute()
tray.Finish()
