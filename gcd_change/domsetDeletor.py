#!/usr/bin/env python
import numpy as np
import icecube
from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
import os

def changeDFrame(dframe):
	outframe  = icetray.I3Frame(icetray.I3Frame.DetectorStatus)
	for ikey in dframe.keys():
		if ikey != "I3DetectorStatus":
			outframe[ikey] = frame[ikey]
		else:
			outframe[ikey] = icecube.dataclasses.I3DetectorStatus()
			outframe[ikey].daq_configuration_name = frame[ikey].daq_configuration_name
			outframe[ikey].dom_status = frame[ikey].dom_status
			outframe[ikey].end_time = frame[ikey].end_time
			outframe[ikey].start_time = frame[ikey].start_time
			status  = frame[ikey].trigger_status
			status_ = outframe[ikey].trigger_status
	
			for key, trigger in status:
				if key.config_id == 102:
					status_[key] = icecube.dataclasses.I3TriggerStatus()
					status_[key].readout_settings = trigger.readout_settings
					status_[key].identity = trigger.identity
					status_[key].trigger_name = trigger.trigger_name
					status_[key].trigger_settings["threshold"] = trigger.trigger_settings["threshold"]
					status_[key].trigger_settings["timeWindow"] = trigger.trigger_settings["timeWindow"]
			else:
				status_[key] = trigger
	return outframe

oldGCD = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"

i3file = dataio.I3File("/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305ModifiedNoDomSet.i3.gz", "w")

for frame in dataio.I3File(oldGCD,'r'):
  if frame.Stop == icetray.I3Frame.DetectorStatus:
    outframe=changeDFrame(frame)
    i3file.push(outframe)
  else:
    i3file.push(frame)
i3file.close()