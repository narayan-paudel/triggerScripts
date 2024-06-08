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
# parser.add_argument('--GCD',"-g",type=str,
#   default="/home/enpaudel/dataExp/run2023/GCD/PFGCD_Run00138807_Subrun00000000.i3.gz",help="GCD file")
args = parser.parse_args()

import numpy as np
import math
from astropy.time import Time

# GCD = args.GCD
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
outputDir = "/data/sim/IceTop/2023/generated/untriggered/run2023/forbushL2Window/"


refTime = dataclasses.I3Time(2024)
print("refTime",refTime)
t_ref_start = Time("2024-03-24 18:00:00").mjd
t_ref_end = Time("2024-03-25 03:00:00").mjd

# t_ref_start = Time("2024-03-24 23:14:17").mjd
# t_ref_end = Time("2024-03-24 23:14:18").mjd

# t_ref = Time("2024-03-24 20:00:00")
# t_ref = t_ref.to_value('mjd', 'long')
# t_refi3 = dataclasses.I3Time(t_ref)
# print(t_ref,dataclasses.I3Time().set_mod_julian_time_double(t_ref),t_refi3, t_refi3.mod_julian_day)
# print(dataclasses.I3Time(t_ref))

class TriggerRate(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.HG7Triggers = 0
    self.HLCTriggers = 0

  def DAQ(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    if frame.Has("I3EventHeader"):
      self.runID = frame["I3EventHeader"].run_id
      # if frame.Has("DSTTriggers"):
        # for trigger in frame['DSTTriggers'].unpack(frame['I3DetectorStatus']):
        #     # if trigger.key.config_id == 30043 and trigger.fired:
        #     if trigger.key.config_id == 102 and trigger.fired:
        #         self.HLCTriggers += 1
        #     if trigger.key.config_id == 30043 and trigger.fired:
        #         self.HG7Triggers += 1


def selectTime(frame):
  inWindow = False
  if frame.Has("I3EventHeader"):
    runID = frame["I3EventHeader"].run_id
    evh = frame["I3EventHeader"]
    start_time = evh.start_time.mod_julian_day+evh.start_time.mod_julian_sec/(24.*3600.)+(evh.start_time.mod_julian_nano_sec*math.pow(10,-9))/(24.*3600.)
    end_time = evh.end_time.mod_julian_day+evh.end_time.mod_julian_sec/(24.*3600.)+(evh.end_time.mod_julian_nano_sec*math.pow(10,-9))/(24.*3600.)
    if start_time >= t_ref_start and start_time <= t_ref_end:
      inWindow = True
      # print("dateTime",start_time.date_time)
      # print("I3Time",dataclasses.I3Time(start_time))
      # time=evh.start_time.mod_julian_day+evh.start_time.mod_julian_sec/(24.*3600.)+(evh.start_time.mod_julian_nano_sec*math.pow(10,-9))/(24.*3600.)
      # print(dataclasses.I3Time().set_mod_julian_time_double(start_time))
      # print(end_time - start_time)
      # print("daq",evh.start_time.utc_daq_time)
      # print("year",evh.start_time.utc_year)
      # print(start_time.mod_julian_day)
      # print(evh.start_time.mod_julian_day+evh.start_time.mod_julian_sec/(24.*3600.)+(evh.start_time.mod_julian_nano_sec*math.pow(10,-9))/(24.*3600.))
  return inWindow





tray = I3Tray()
tray.AddModule("I3Reader","reader",
              # FilenameList=[GCD]+inputList,
              FilenameList=inputList,
              # filenameList=inputList[0],
              # filename=GCD,
              )
# tray.AddModule(TriggerRate, "rT",
#             # GCD=GCD,
#             # Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
#             )
tray.AddModule(selectTime, "ITTrig",
            # GCD=GCD,
            Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
            )

tray.AddModule("I3Writer","i3writer",
            filename=str(outputDir)+filename.split("/")[-1].split(".")[0] + "Win.i3.gz",
            streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            # DropOrphanStreams=[icetray.I3Frame.Physics]
            )

tray.Execute()
tray.Finish()