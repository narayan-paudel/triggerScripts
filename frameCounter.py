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

# failedFile = /home/enpaudel/icecube/triggerStudy/simFiles/dataSetGen/FeDAT001977Gen.i3.bz2

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
		print(args.input)
		if int(self.eventNo) != 100:
			print("Some events are lost, recheck scripts",self.eventNo)
			with open("../lessThan100.txt","a") as f:
				f.write(str(*args.input))
				f.write(" "+str(self.eventNo)+"\n")

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