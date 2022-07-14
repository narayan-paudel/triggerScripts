#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

from scipy.spatial import ConvexHull

from icecube import dataclasses, icetray, dataio, radcube, rock_bottom
from I3Tray import *

import numpy as np

import os
ABS_PATH_HERE=str(os.path.dirname(os.path.realpath(__file__)))
outputDir = "/home/enpaudel/icecube/triggerStudy/simFiles/"

# args = parser.parse_args()

GCDFile_Kath = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"
GCDFile_Alan = "/home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz"
framenames=["I3AntennaGeometry", "I3Geometry", "I3ScintGeometry", "I3IceActGeometry"]

f = dataio.I3File(GCDFile_Kath)
gframe = f.pop_frame()
while f.more() and not "I3Geometry" in gframe:
    gframe = f.pop_frame()

def snowCopy(gframe):
	I3Geo = gframe["I3Geometry"].stationgeo
	for key in I3Geo.keys():
		station = I3Geo[key]
		for tank in station:
			pos = tank.position
			snowHeight = tank.snowheight
			# print(snowHeight)
	# print(I3Geo)

# snowCopy(gframe)


def removeITVolume(frame):
	if frame.Has("I3DetectorStatus"):
		detStatus = frame["I3DetectorStatus"]
	for key, trigger in detStatus.trigger_status.items():
		print (key.config_id, key.source.name, key.type.name, trigger.trigger_settings.get('threshold', 'N/A'),
		 [(source.name, config.readout_time_minus, config.readout_time_plus) for source, config in trigger.readout_settings.items()])
	trigStatus = detStatus.trigger_status
	trigStatNew = []
	for key,trigger in trigStatus:
		if key.config_id == 21002:
			del trigStatus[key]
	print("trigstatus",trigStatNew)
	# trigStatusNew = dataclasses.I3DetectorStatus.trigger_status(trigStatNew)
	del frame["I3DetectorStatus"]
	detStatus.trigger_status = trigStatus
	frame["I3DetectorStatus"] = detStatus
	# trigStatusNew = dataclasses.I3TriggerStatus()
	# selectTriggerStatus = [(key,ts) for key,ts in detStatus.trigger_status.items() if key.config_id != 21002]
	# print("selected triggers",selectTriggerStatus)
	# print(detStatus.trigger_status.items())
	# print(detStatus.trigger_status)
	for key, trigger in detStatus.trigger_status.items():
		print (key.config_id, key.source.name, key.type.name, trigger.trigger_settings.get('threshold', 'N/A'),
		 [(source.name, config.readout_time_minus, config.readout_time_plus) for source, config in trigger.readout_settings.items()])
def removeDOMSet_(trigStatus):
	# newTrigSetting = dataclasses.I3TriggerStatus().trigger_settings["domSet"]
	for key,trigger in trigStatus:
		if key.config_id == 102:
			smtKey = key
	trigger = trigStatus[smtKey].trigger_settings["domSet"]
	print("trigger domset",trigger)
	del trigStatus[smtKey].trigger_settings["domSet"]
	return trigStatus

def removeDOMSet(frame):
	if frame.Has("I3DetectorStatus"):
		detStatus = frame["I3DetectorStatus"]
	for key, trigger in detStatus.trigger_status.items():
		print (key.config_id, key.source.name, key.type.name, trigger.trigger_settings.get('threshold', 'N/A'),
		 [(source.name, config.readout_time_minus, config.readout_time_plus) for source, config in trigger.readout_settings.items()])
	trigStatus = detStatus.trigger_status
	trigStatNew = []
	trigStatus = removeDOMSet_(trigStatus)
	# for key,trigger in trigStatus:
	# 	if key.config_id == 102:
	# 		# print("SMT trigger",key,trigger,trigger.trigger_settings["domSet"])
	# 		print("SMT trigger domset",key,trigger.trigger_settings["domSet"])
	# 		del trigger.trigger_settings["domSet"]
	# 		# print("SMT trigger domset",key,trigger.trigger_settings["domSet"])
	# for key,trigger in trigStatus:
	# 	if key.config_id == 102:
	# 		# print("SMT trigger",key,trigger,trigger.trigger_settings["domSet"])
	# 		print("SMT trigger domset",key,trigger.trigger_settings["domSet"])
	# 		del trigger.trigger_settings["domSet"]
	print("trigstatus",trigStatNew)
	# trigStatusNew = dataclasses.I3DetectorStatus.trigger_status(trigStatNew)
	del frame["I3DetectorStatus"]
	detStatus.trigger_status = trigStatus
	frame["I3DetectorStatus"] = detStatus
	# trigStatusNew = dataclasses.I3TriggerStatus()
	# selectTriggerStatus = [(key,ts) for key,ts in detStatus.trigger_status.items() if key.config_id != 21002]
	# print("selected triggers",selectTriggerStatus)
	# print(detStatus.trigger_status.items())
	# print(detStatus.trigger_status)
	for key, trigger in detStatus.trigger_status.items():
		print (key.config_id, key.source.name, key.type.name, trigger.trigger_settings.get('threshold', 'N/A'),
		 [(source.name, config.readout_time_minus, config.readout_time_plus) for source, config in trigger.readout_settings.items()])

tray = I3Tray()
tray.AddModule("I3Reader","reader",
	           filenameList=[GCDFile_Kath],
	           # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )
# tray.AddModule(removeITVolume,"ITVolume",
# 				streams=[icetray.I3Frame.DetectorStatus])

tray.AddModule(removeDOMSet,"domset",
				streams=[icetray.I3Frame.DetectorStatus])

tray.AddModule("I3Writer","i3writer",
	          # filename=str(outputDir)+"/hadronTimeTest/"+str(fileName)+"HadronNonTrigImProperTEvts.i3.gz",
	          # filename=str(outputDir)+"/dataSetUnique/"+str(fileName)+"Unique.i3.gz",
	          filename=str(outputDir)+"/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz",
	          # streams=[icetray.I3Frame.TrayInfo,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	          streams=[icetray.I3Frame.TrayInfo,icetray.I3Frame.Geometry,icetray.I3Frame.Calibration,icetray.I3Frame.DetectorStatus,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	          # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	          # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )
tray.Execute()
tray.Finish()