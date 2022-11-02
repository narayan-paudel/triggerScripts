#!/usr/bin/env python


import os
import sys
from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
import argparse
import math


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/DAT059871GenDet.i3.bz2", help='Input files after running detector.py.')
args = parser.parse_args()


class FrameCounter(icetray.I3Module):
	def __init__(self,ctx):
		icetray.I3Module.__init__(self,ctx)
	def Configure(self):
		self.eventNo = 0

	def DAQ(self,frame):
		if frame.Stop == icetray.I3Frame.DAQ:
			self.eventNo += 1

	def Finish(self):
		# print("total no of events in the frame is ", self.eventNo)
		if int(self.eventNo) != 100:
			print("Some events are lost, recheck scripts",self.eventNo)
			print(args.input)

tray = I3Tray()
tray.AddModule("I3Reader","reader",
	           filenameList=args.input,
	           # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )

tray.AddModule(FrameCounter,"frameCount",
	           # Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           	)
tray.Execute()
tray.Finish()