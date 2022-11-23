#!/usr/bin/env python3

from I3Tray import I3Tray,I3Units
from icecube import dataclasses,icetray,dataio


GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"

f = dataio.I3File(GCD)
dframe = f.pop_frame()
while f.more() and not "I3DetectorStatus" in dframe:
    dframe = f.pop_frame()

def addTankTrigger(frame):
  outframe = icetray.I3Frame(icetray.I3Frame.DetectorStatus)
  for ikey in frame.keys():
    if ikey != "I3DetectorStatus":
      outframe[ikey] = frame[ikey]
    else:
      outframe[ikey] = dataclasses.I3DetectorStatus()
      outframe[ikey].daq_configuration_name = frame[ikey].daq_configuration_name
      outframe[ikey].dom_status = frame[ikey].dom_status
      outframe[ikey].end_time = frame[ikey].end_time
      outframe[ikey].start_time = frame[ikey].start_time
      status  = frame[ikey].trigger_status
      status_ = outframe[ikey].trigger_status
      for key, trigger in status:
        if key.config_id == 102:
          keyTank = dataclasses.TriggerKey(key)
          keyTank.config_id = 108
          triggerTank = dataclasses.I3TriggerStatus(trigger)
      triggerTank.trigger_settings["threshold"] = str(6)
      triggerTank.trigger_settings["timeWindow"] = str(6000)
      status_[keyTank] = triggerTank
      for key, trigger in status:
        status_[key] = trigger
      # print(triggerTank.trigger_settings["threshold"])
      # print(triggerTank.trigger_settings["timeWindow"])
      # # print("I3TriggerReadoutConfig",trigger.I3TriggerReadoutConfig)
      # # print("readout settings",trigger.readout_settings)
      # # print("subdetector",trigger.Subdetector)
      # # print("trigger_name",trigger.trigger_name)
      # # print("trigger_settings",trigger.trigger_settings)
      # for keys, value in trigger.trigger_settings:
      #   print("keys",keys,type(keys),value,type(value))
      # print("trigger threshold",trigger.trigger_settings['threshold'])
      # # trigger.trigger_settings['threshold'] = 5
      # print("trigger timeWindow",trigger.trigger_settings['timeWindow'])
      # # print("key",key,keyTank)
      # status_[keyTank].readout_settings = triggerTank.readout_settings
      # status_[keyTank].identity = triggerTank.identity
      # status_[keyTank].trigger_name = triggerTank.trigger_name
      # status_[keyTank].trigger_settings["timeWindow"] = 5000
  return outframe

def checkTankTrigger(frame):
  detStatus = frame["I3DetectorStatus"]
  for key,trigger in detStatus.trigger_status:
    if key.config_id == 108 or key.config_id == 102:
      print("key",key)
      print("threshold",trigger.trigger_settings["threshold"])
      print("timeWindow",trigger.trigger_settings["timeWindow"])
      print("I3TriggerReadoutConfig",trigger.I3TriggerReadoutConfig)
      print("readout settings",trigger.readout_settings)
      print("subdetector",trigger.Subdetector)
      print("trigger_name",trigger.trigger_name)
      print("trigger_settings",trigger.trigger_settings)
      print("trigStatus",trigger,trigger.trigger_settings.get('threshold', 'N/A'))
      print("trigStatus",trigger,trigger.trigger_settings.get('timeWindow', 'N/A'))
      print("trigStatus",trigger,trigger.trigger_settings.get('domSet', 'N/A'))




GCDTankTrig = dataio.I3File("/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz", "w")
for frame in dataio.I3File(GCD,'r'):
  if frame.Stop == icetray.I3Frame.DetectorStatus:
    outframe=addTankTrigger(frame)
    GCDTankTrig.push(outframe)
  else:
    GCDTankTrig.push(frame)
GCDTankTrig.close()

for frame in dataio.I3File("/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz", "r"):
  if frame.Stop == icetray.I3Frame.DetectorStatus:
    checkTankTrigger(frame)
