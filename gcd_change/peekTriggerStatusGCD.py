#!/usr/bin/env python3

from I3Tray import I3Tray,I3Units
from icecube import dataclasses,icetray,dataio


GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"

f = dataio.I3File(GCD)
dframe = f.pop_frame()
while f.more() and not "I3DetectorStatus" in dframe:
    dframe = f.pop_frame()

tankConfigIDMap = {103:(6,5000),104:(7,5000),105:(8,5000),106:(9,5000),107:(10,5000),
113:(6,4000),114:(7,4000),115:(8,4000),116:(9,4000),117:(10,4000),
123:(6,3000),124:(7,3000),125:(8,3000),126:(9,3000),127:(10,3000),
133:(6,2000),134:(7,2000),135:(8,2000),136:(9,2000),137:(10,2000)} #map of config ID and no of tanks

def updateConfigID(smtkey,configID):
  """
  updates config ID of given smtkey
  """
  keyTank = dataclasses.TriggerKey(smtkey)
  keyTank.config_id = configID
  return keyTank


def addTrigger(status,tankConfigIDMap,keyTank,triggerTank):
  for iConfigID in tankConfigIDMap.keys():
    keyTank.config_id = iConfigID
    triggerTank.trigger_settings["threshold"] = str(tankConfigIDMap[iConfigID][0])
    triggerTank.trigger_settings["timeWindow"] = str(tankConfigIDMap[iConfigID][1])
    status[keyTank] = triggerTank
  return status

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
        status_[key] = trigger
        if key.config_id == 102:
          smtkey = dataclasses.TriggerKey(key)
          smtTrigger = dataclasses.I3TriggerStatus(trigger)
      status_ = addTrigger(status_,tankConfigIDMap,smtkey,smtTrigger)
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
    if key.config_id in [103,104,105,106,107,113,114,115,116,117,123,124,125,126,127,133,134,135,136,137]:
      print("key",key)
      print("threshold",trigger.trigger_settings["threshold"])
      print("timeWindow",trigger.trigger_settings["timeWindow"])
      print("I3TriggerReadoutConfig",trigger.I3TriggerReadoutConfig)
      print("readout settings",trigger.readout_settings)
      print("subdetector",trigger.Subdetector)
      print("trigger_name",trigger.trigger_name, type(trigger.trigger_name))
      print("trigger_settings",trigger.trigger_settings)
      # print("trigStatus",trigger,trigger.trigger_settings.get('threshold', 'N/A'))
      # print("trigStatus",trigger,trigger.trigger_settings.get('timeWindow', 'N/A'))
      # print("trigStatus",trigger,trigger.trigger_settings.get('domSet', 'N/A'))

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
