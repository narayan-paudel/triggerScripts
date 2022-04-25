#!/usr/bin/env python
import numpy as np
from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
import os

gcd_file2020 = "/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz"
gcd_file2016 = "/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2016.57531_V0_OctSnow.i3.gz"
gcd_file2012 = "/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz"
gcd_file2011 = "/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2011.55697_V3.i3.gz"
gcd_fileAlan = "/home/enpaudel/icecube/triggerStudy/triggerScripts/../simFiles/GCD-Survey-AntITScint_2020.02.24.i3.gz"
gcd_fileKath = "/home/enpaudel/icecube/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"
gcd_fileMy = "/home/enpaudel/icecube/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"
gcd_fileAgni = "/home/enpaudel/icecube/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305ModifiedNoDomSet.i3.gz"

# gcd2020_i3 = dataio.I3File(gcd_file2020)
# gcd2016_i3 = dataio.I3File(gcd_file2016,'r')

def checkDOMSet(gcd):
	f = dataio.I3File(gcd)	
	dframe = f.pop_frame()
	while f.more() and not "I3DetectorStatus" in dframe:
	    dframe = f.pop_frame()
	trigStatus = dframe["I3DetectorStatus"].trigger_status
	for key,trigger in trigStatus:
		if key.config_id == 102:
			smtKey = key
			print(trigger.trigger_settings["domSet"])
	dmset = trigStatus[smtKey].trigger_settings["domSet"]
	print("dmset",dmset)


checkDOMSet(gcd_file2020)
checkDOMSet(gcd_fileKath)
# checkDOMSet(gcd_fileMy)
# checkDOMSet(gcd_fileAgni)
# checkDOMSet(gcd_file2016)
# checkDOMSet(gcd_file2012)
# checkDOMSet(gcd_file2011)
checkDOMSet(gcd_fileAlan)
