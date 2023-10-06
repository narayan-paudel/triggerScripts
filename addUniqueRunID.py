#!/usr/bin/env python

import os
import sys
import subprocess

from I3Tray import I3Tray, I3Units
from icecube.simprod import segments
from icecube.simprod.util import simprodtray, arguments
from icecube.simprod.util import ReadI3Summary, WriteI3Summary
from icecube import icetray, dataclasses, dataio
from icecube.simprod.util.simprodtray import RunI3Tray
import argparse
import icecube.icetray
import math
from icecube import phys_services, sim_services
from icecube import tableio, hdfwriter
from icecube.simprod.util import PrintContext
from icecube import topeventcleaning, tpx


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input',"-i", type=str,
 default="/home/enpaudel/icecube/triggerStudy/simFiles/ITGenTest/FeDAT000001GenDet.i3.bz2", help='Input files after running detector.py.')
parser.add_argument('--output',"-o", type=str,
 default="/home/enpaudel/icecube/triggerStudy/simFiles/ITGenTest/FeDAT000001GenDetUnique.i3.bz2", help='Output files after running detector.py.')
parser.add_argument('--GCD',"-g", type=str,
 default="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz", help='Output files after running detector.py.')
args = parser.parse_args()

def addUniqieEvtID(frame):
	newHead = dataclasses.I3EventHeader()
	if frame.Has("MCPrimary"):
		particle_type = str(frame["MCPrimary"].type)
		if "PP" in particle_type:
			particleID = 1
		elif "He" in particle_type:
			particleID = 2
		elif "O" in particle_type:
			particleID = 3
		elif "Fe" in particle_type:
			particleID = 4
		# print("particle type",particle_type)
	if frame.Has("I3EventHeader"):
		newHead.start_time = frame["I3EventHeader"].start_time
		newHead.end_time = frame["I3EventHeader"].end_time
		newHead.sub_event_stream = frame["I3EventHeader"].sub_event_stream
		newHead.sub_event_id = frame["I3EventHeader"].sub_event_id
		newHead.sub_run_id = frame["I3EventHeader"].sub_run_id
		runID = frame["I3EventHeader"].run_id
		evtID = frame["I3EventHeader"].event_id
		newEvtID = runID * 1000 + particleID * 100 + evtID
		newHead.run_id = runID
		newHead.event_id = newEvtID
		del frame["I3EventHeader"]
		frame.Put("I3EventHeader",newHead)
		# print("new event header",frame["I3EventHeader"].event_id)


tray = I3Tray()
tray.AddModule("I3Reader","reader",
	           filenameList=[args.GCD]+[args.input],
	          )
tray.AddModule(addUniqieEvtID,"uniqEvtID",
				streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

tray.AddModule("I3Writer","i3writer",
	          filename=args.output,
	          streams=[icetray.I3Frame.TrayInfo,icetray.I3Frame.DAQ,icetray.I3Frame.Physics,icetray.I3Frame.Simulation],
	          )
tray.Execute()
tray.Finish()