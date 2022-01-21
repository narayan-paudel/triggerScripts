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
args = parser.parse_args()

fileName = args.input[0].split("/")[-1]
fileName = fileName.split(".")[0]
print(fileName)


class TriggerCheck(icetray.I3Module):
	def __init__(self,ctx):
		icetray.I3Module.__init__(self,ctx)
	def Configure(self):
		self.SMTTriggered = 0
		self.SMTUntriggered = 0
		self.GlobalTriggered = 0
		self.GlobalUntriggered = 0
		self.totalEvents = 0


	def DAQ(self,frame):
		"""what if there is only one trigger in trigger hierarchy
		"""
		# trigger = frame["I3Triggers"]
		triggerHierarchy = frame["QTriggerHierarchy"]
		# print(frame.keys())
		self.totalEvents += 1
		print("total no of events",self.totalEvents)
		print("triggerHierarchy length",len(triggerHierarchy))
		if len(triggerHierarchy) == 0:
			frame["ITSMTTriggered"] = dataclasses.I3Double(0)
			frame["ITGlobalTriggered"] = dataclasses.I3Double(0)
			self.GlobalUntriggered+=1
			self.SMTUntriggered+=1
		elif len(triggerHierarchy) == 1:
			for t in triggerHierarchy:
				if (t.key.config_id == 102):
					frame["ITGlobalTriggered"] = dataclasses.I3Double(0)
					self.GlobalUntriggered+=1
					if t.fired:
						frame["ITSMTTriggered"] = dataclasses.I3Double(1)
						self.SMTTriggered+=1
					elif not t.fired:
						frame["ITSMTTriggered"] = dataclasses.I3Double(0)
						self.SMTUntriggered+=1
				elif (t.key.config_id == None):					
					frame["ITSMTTriggered"] = dataclasses.I3Double(0)
					self.SMTUntriggered+=1
					if t.fired:
						frame["ITGlobalTriggered"] = dataclasses.I3Double(1)
						self.GlobalTriggered+=1
					elif not t.fired:
						frame["ITGlobalTriggered"] = dataclasses.I3Double(0)
						self.GlobalUnTriggered+=1
		else:
			# print([t.key.config_id for t in triggerHierarchy])			
			for t in triggerHierarchy:
				# print("tk",t.key)
				if (t.key.config_id == 102 and t.fired):
					frame["ITSMTTriggered"] = dataclasses.I3Double(1)
					self.SMTTriggered+=1
				elif (t.key.config_id == 102 and not t.fired):
					frame["ITSMTTriggered"] = dataclasses.I3Double(0)
					self.SMTUntriggered+=1
				if (t.key.config_id == None and t.fired):
					frame["ITGlobalTriggered"] = dataclasses.I3Double(1)
					self.GlobalTriggered+=1
				elif (t.key.config_id == None and not t.fired):
					frame["ITGlobalTriggered"] = dataclasses.I3Double(0)
					self.GlobalUntriggered+=1
		self.PushFrame(frame)

	def Finish(self):
		print("total number of events",self.totalEvents)
		print("number of SMT triggered events",self.SMTTriggered )
		print("number of SMT untriggered events",self.SMTUntriggered )
		print("number of Global triggered events",self.GlobalTriggered)
		print("number of Gloabl untriggered events",self.GlobalUntriggered)
		# hdftable = hdfwriter.I3HDFTableService(str(outputDir)+str(fileName)+"TrigCount.hdf5")




def TriggerCheck2(frame):
	trigger = frame["I3Triggers"]
	triggerHierarchy = frame["QTriggerHierarchy"]
	# print("trigger",trigger)
	# print("hierarchy",triggerHierarchy)
	# triggers = dataclasses.I3TriggerHierarchy.from_frame(frame, key)
	# a=0.0
	for t in triggerHierarchy:
		# print("tk",t.key)
		if (t.key.config_id == 102 and t.fired):
			frame["ITSMTTriggered"] = dataclasses.I3Double(1)
		elif (t.key.config_id == 102 and not t.fired):
			frame["ITSMTTriggered"] = dataclasses.I3Double(0)
		if (t.key.config_id == None and t.fired):
			frame["ITGlobalTriggered"] = dataclasses.I3Double(1)
		elif (t.key.config_id == None and not t.fired):
			frame["ITGlobalTriggered"] = dataclasses.I3Double(0)


def AddTotalCharge(frame):
	'''calculates total SLC charge in SLC tank pulses'''
	# print("Anything")
	if (frame.Has('OfflineIceTopSLCTankPulses')):
		# print("frame has OfflineIceTopSLCTankPulses")
		slc = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'OfflineIceTopSLCTankPulses')
		ITTotalChargeSLC = sum([sum([c.charge for c in slc[i]]) for i in slc.keys()])
		frame["ITTotalChargeSLC"] = dataclasses.I3Double(ITTotalChargeSLC)

def AddTotalCharge(frame,keys):
	'''calculates total SLC or HLC charges in tank pulses
	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	'''
	# print("Anything")
	for key in keys:		
		if (frame.Has(str(key))):
			# print("frame has", str(key))
			lc = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,str(key))
			ITTotalCharge = sum([sum([c.charge for c in lc[i]]) for i in lc.keys()])
			frame[str(key)+"TotalCharge"] = dataclasses.I3Double(ITTotalCharge)

def AddTotalTankHit(frame,pulseseriesList):
	'''calculates total SLC or HLC charges in tank pulses
	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	'''
	# print("Adding tank Hit")
	for pulseseries in pulseseriesList:
		NCh = 0
		psm = frame[pulseseries]
		if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
		for om,pulses in psm:
			for pulse in pulses:
				NCh +=1
				break
		channels = [omkey for omkey,ps in psm if len(ps)>0]
		if NCh == 0:
			print("psm",psm)
		frame[str(pulseseries)+"TotalTankHit"] = dataclasses.I3Double(NCh)
		print("No. of Hits in ",pulseseries,NCh,len(channels),frame[str(pulseseries)+"TotalTankHit"])



tray = I3Tray()
tray.AddModule("I3Reader","reader",
	           filenameList=[GCD]+args.input,
	           # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )
tray.AddModule('I3HLCTankPulseMerger',
	           # filenameList=[GCD]+args.input,
	           InputVEMPulses = 'OfflineIceTopSLCVEMPulses',
	           OutputTankPulses = 'OfflineIceTopSLCTankPulses',
	           ExcludedTanks  = 'ExcludedSLCTanks')

tray.AddModule(TriggerCheck, "trigChk",
				# Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
				)
# 
# tray.AddModule(AddTotalCharge,"addCharge",
# 	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
# 	           	)

tray.AddModule(AddTotalCharge,"addCharge",
	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           keys=['OfflineIceTopSLCTankPulses','OfflineIceTopHLCTankPulses']
	           	)
tray.AddModule(AddTotalTankHit,"addHit",
	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           pulseseriesList=['OfflineIceTopSLCTankPulses','OfflineIceTopHLCTankPulses']
	           # pulseseriesList=['OfflineIceTopHLCTankPulses']
	           	)

tray.AddModule("I3Writer","i3writer",
	          filename=str(outputDir)+"trigChkTest/"+str(fileName)+"TrigCount.i3.gz",
	          # streams=[icetray.I3Frame.DAQ],
	          streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	          # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )

# tray.AddModule("I3NullSplitter","nullsplitter")

# tray.Add(hdfwriter.I3HDFWriter, 'hdf',
#     Output=str(outputDir)+"/trigChkTest/"+str(fileName)+"TrigCount.hdf5",
#     CompressionLevel=9,
#     # SubEventStreams=['IceTopSplit'],
#     SubEventStreams=['NullSplit'],
#     # SubEventStreams=["nullsplitter",'IceTopSplit',"nullsplitter",'NullSplit',]
#     # SubEventStreams=["ice_top"],
#     # Streams=[icetray.I3Frame.DAQ],
#     keys = ["MCPrimary","ITSMTTriggered","ITGlobalTriggered","OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalTankHit","OfflineIceTopSLCTankPulsesTotalTankHit"],
#     )

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