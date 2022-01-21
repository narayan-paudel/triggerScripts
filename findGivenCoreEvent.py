#!/usr/bin/env python


import os
import sys
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


CORSIKA_ID = "DAT059871"
outputDir = "/home/enpaudel/icecube/triggerStudy/simFiles/"

GCD = "/home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz"

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/DAT059871GenDet.i3.bz2", help='Input files after running detector.py.')
parser.add_argument('--xCore',"-x",type=float,default=0,help='x cordinate of required core')
parser.add_argument('--yCore',"-y",type=float,default=0,help='x cordinate of required core')
args = parser.parse_args()

fileName = args.input[0].split("/")[-1]
fileName = fileName.split(".")[0]
print(fileName)


nFrame = 0


def findCore(frame):
	'''calculates total SLC charge in SLC tank pulses'''
	# print("Anything")
	global nFrame
	if (frame.Has("MCPrimary")):
		if abs(frame["MCPrimary"].pos.x - args.xCore) < 0.5 and abs(frame["MCPrimary"].pos.y - args.yCore) < 0.5 :
			print("yell there was the problem",nFrame,frame["MCPrimary"].pos.x,frame["MCPrimary"].pos.y)
	nFrame += 1







tray = I3Tray()
tray.AddModule("I3Reader","reader",
	           filenameList=[GCD]+args.input,
	           # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )

tray.AddModule(findCore,"coreFinder",
	streams=[icetray.I3Frame.DAQ])


tray.AddModule("I3Writer","i3writer",
	          filename=str(outputDir)+"timingDebug/timingDebugCore.i3.gz",
	          # streams=[icetray.I3Frame.DAQ],
	          streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	          # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )

tray.Execute()
tray.Finish()



#################################################333
###################################################
# hdftable = hdfwriter.I3HDFTableService(str(outputDir)+"/trigChk/"+str(fileName)+"TrigCount.hdf5")
# tray.AddModule(tableio.I3TableWriter,'hdf1',
# 	         tableservice = hdftable,
# 	         # SubEventStreams = ['in_ice'],
# 	         keys = ['I3EventHeader',"I3Triggers","ITSMTTriggered","ITGlobalTriggered"],
# 	         # types = [dataclasses.I3Particle] #inplace of keys
# 	        )