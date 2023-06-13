#!/usr/bin/env python3

from I3Tray import I3Tray,I3Units
from icecube import dataclasses,icetray,dataio


GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"

f = dataio.I3File(GCD)
dframe = f.pop_frame()
while f.more() and not "I3DetectorStatus" in dframe:
    dframe = f.pop_frame()
    
# tankConfigIDMap = {103:(6,5000),104:(7,5000),105:(8,5000),106:(9,5000),107:(10,5000),
# 113:(6,4000),114:(7,4000),115:(8,4000),116:(9,4000),117:(10,4000),
# 123:(6,3000),124:(7,3000),125:(8,3000),126:(9,3000),127:(10,3000),
# 133:(6,2000),134:(7,2000),135:(8,2000),136:(9,2000),137:(10,2000)} # initial map of config ID and no of tanks

tankConfigIDMap = {103:(4,5000,"HLC"),104:(2,5000,"HLC"),
161:(6,1000,"both"),162:(6,2000,"both"),163:(6,3000,"both"),164:(6,4000,"both"),165:(6,5000,"both"),
171:(7,1000,"both"),172:(7,2000,"both"),173:(7,3000,"both"),174:(7,4000,"both"),175:(7,5000,"both"),
181:(8,1000,"both"),182:(8,2000,"both"),183:(8,3000,"both"),184:(8,4000,"both"),185:(8,5000,"both"),
191:(9,1000,"both"),192:(9,2000,"both"),193:(9,3000,"both"),194:(9,4000,"both"),195:(9,5000,"both"),
1101:(10,1000,"both"),1102:(10,2000,"both"),1103:(10,3000,"both"),1104:(10,4000,"both"),1105:(10,5000,"both"),
261:(6,1000,"both_hg"),262:(6,2000,"both_hg"),263:(6,3000,"both_hg"),264:(6,4000,"both_hg"),265:(6,5000,"both_hg"),
271:(7,1000,"both_hg"),272:(7,2000,"both_hg"),273:(7,3000,"both_hg"),274:(7,4000,"both_hg"),275:(7,5000,"both_hg"),
281:(8,1000,"both_hg"),282:(8,2000,"both_hg"),283:(8,3000,"both_hg"),284:(8,4000,"both_hg"),285:(8,5000,"both_hg"),
291:(9,1000,"both_hg"),292:(9,2000,"both_hg"),293:(9,3000,"both_hg"),294:(9,4000,"both_hg"),295:(9,5000,"both_hg"),
2101:(10,1000,"both_hg"),2102:(10,2000,"both_hg"),2103:(10,3000,"both_hg"),2104:(10,4000,"both_hg"),2105:(10,5000,"both_hg"),
321:(2,1000,"SLC"),322:(2,2000,"SLC"),323:(2,3000,"SLC"),324:(2,4000,"SLC"),325:(2,5000,"SLC"),
341:(4,1000,"SLC"),342:(4,2000,"SLC"),343:(4,3000,"SLC"),344:(4,4000,"SLC"),345:(4,5000,"SLC"),
361:(6,1000,"SLC"),362:(6,2000,"SLC"),363:(6,3000,"SLC"),364:(6,4000,"SLC"),365:(6,5000,"SLC"),
371:(7,1000,"SLC"),372:(7,2000,"SLC"),373:(7,3000,"SLC"),374:(7,4000,"SLC"),375:(7,5000,"SLC"),
381:(8,1000,"SLC"),382:(8,2000,"SLC"),383:(8,3000,"SLC"),384:(8,4000,"SLC"),385:(8,5000,"SLC"),
391:(9,1000,"SLC"),392:(9,2000,"SLC"),393:(9,3000,"SLC"),394:(9,4000,"SLC"),395:(9,5000,"SLC"),
3101:(10,1000,"SLC"),3102:(10,2000,"SLC"),3103:(10,3000,"SLC"),3104:(10,4000,"SLC"),3105:(10,5000,"SLC")
}

# tankConfigIDMap = {103:(4,5000,"HLC"),104:(2,5000,"HLC"),
# }

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
    keyTank.source = dataclasses.ICE_TOP
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
    # if key.config_id in [103,104,105,106,107,113,114,115,116,117,123,124,125,126,127,133,134,135,136,137]:
    if key.config_id in tankConfigIDMap.keys():
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

# GCDTankTrig = dataio.I3File("/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz", "w")
GCDTankTrig = dataio.I3File("/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/../GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrigAll.i3.gz", "w")
for frame in dataio.I3File(GCD,'r'):
  if frame.Stop == icetray.I3Frame.DetectorStatus:
    outframe=addTankTrigger(frame)
    GCDTankTrig.push(outframe)
  else:
    GCDTankTrig.push(frame)
GCDTankTrig.close()

for frame in dataio.I3File("/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz", "r"):
# for frame in dataio.I3File("/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/../GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz", "r"):
  if frame.Stop == icetray.I3Frame.DetectorStatus:
    checkTankTrigger(frame)
