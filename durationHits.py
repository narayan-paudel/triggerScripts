#!/usr/bin/env python3


import os
import sys
from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube import topeventcleaning, tpx,toprec

import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+',
 default="/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTanks/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz",
 help='Input files after running detector.py.')
args = parser.parse_args()

GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"

exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

class checkDeltaT(icetray.I3Module):
	def __init__(self,ctx):
		icetray.I3Module.__init__(self,ctx)
	def Configure(self):
		self.SMTTriggered = 0
		self.SMTUntriggered = 0
		self.GlobalTriggered = 0
		self.GlobalUntriggered = 0
		self.totalEvents = 0
		self.delta_t_lowEVert = []
		self.delta_t_lowEIncl = []
		self.delta_t_highEVert = []
		self.delta_t_highEIncl = []


	def checkTriggers(self,frame,triggerHierarchy,config_id):
		trigLabel = ConfigIDMap[config_id][2]+"_"+str(ConfigIDMap[config_id][1])
		if len(triggerHierarchy) == 0:
			frame[trigLabel] = dataclasses.I3Double(0)
		else:
			SMTTriggers = [t for t in triggerHierarchy if (t.key.config_id == config_id and t.fired)]
			if len(SMTTriggers) != 0:
				frame[trigLabel] = dataclasses.I3Double(1)
			else:
				frame[trigLabel] = dataclasses.I3Double(0)
		return frame
	def deltaTHit(self,frame,pulseseries):
		psm = frame[pulseseries]
		if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
		times_list = []
		hit_stations = []
		hit_omkeys = []
		for om,pulses in psm:
			hit_stations.append(om[0])
			hit_omkeys.append(om)
			break
		hit_stations = list(set(hit_stations))
		hit_tanks = []
		hit_times = []
		pulsesList = []
		for istation in hit_stations:
			tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
			(iom[0] == istation and iom[0] in [26,39,74] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
			pulsesList += [ipulse for om in tanks for pulse in psm[om] for ipulse in pulse]
		pulseList = sorted(pulsesList, key=lambda pulse: pulse.time)
		hit_times = [pulse.time for pulse in pulseList]
		deltaT = np.array([])
		if len(hit_times) > 0:
			deltaT = np.diff(hit_times)
		return deltaT

	def DAQ(self,frame):
		"""what if there is only one trigger in trigger hierarchy
		"""
		if frame["SMT273"] > 0:
			if np.rad2deg(frame["MCPrimary"].dir.zenith) <= 45 and 10**6.2<=frame["MCPrimary"].energy<10**6.8:
				self.delta_t_lowEVert += self.deltaTHit(frame["IceTopVEMPulses"])
			if np.rad2deg(frame["MCPrimary"].dir.zenith) > 45 and 10**6.2<=frame["MCPrimary"].energy<10**6.8:
				self.delta_t_lowEIncl += self.deltaTHit(frame["IceTopVEMPulses"])
			if np.rad2deg(frame["MCPrimary"].dir.zenith) <= 45 and 10**7.0<=frame["MCPrimary"].energy<10**7.5:
				self.delta_t_highEVert += self.deltaTHit(frame["IceTopVEMPulses"])
			if np.rad2deg(frame["MCPrimary"].dir.zenith) > 45 and 10**7.0<=frame["MCPrimary"].energy<10**7.5:
				self.delta_t_highEIncl += self.deltaTHit(frame["IceTopVEMPulses"])
		# print(frame.keys())
		self.totalEvents += 1
		self.PushFrame(frame)

	def Finish(self):
		print("total number of events",self.totalEvents)


tray = I3Tray()
tray.AddModule("I3Reader","reader",
	filenameList=[GCD]+args.input,
						)
keep_list = ['I3DetectorStatus','I3Geometry','I3Triggers','CleanIceTopRawData','H4aWeight',
  	'HLC6_5000','tank7_3000','I3EventHeader','QTriggerHierarchy','QFilterMask',
	'IceTopRawData', 'MCPrimary','IceTopPulses',
	'OfflineIceTopHLCVEMPulses', 'OfflineIceTopHLCTankPulses','OfflineIceTopSLCTankPulses', 'OfflineIceTopSLCVEMPulses',
	 "SMT183", "SMT263", "SMT273", "SMT283",]

tray.AddModule("Keep","keep unusual time",
 	           keys = keep_list,
 	           # If = True,
 	           streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
 	           )
def Unify(frame, Keys, Output):
  """
  Simple utility to merge RecoPulseSerieses into a single Union.
  """
  extants = [k for k in Keys if k in frame]
  union = dataclasses.I3RecoPulseSeriesMapUnion(frame, extants)
  frame[Output] = union

tray.Add(Unify,"UnionHLCSLC",
  Keys=['OfflineIceTopHLCTankPulses','OfflineIceTopSLCTankPulses'],
  Output='IceTopTankPulses',
  streams=[icetray.I3Frame.DAQ],
  )

tray.Add(Unify,"UnionHLCSLCVEM",
  Keys=['OfflineIceTopHLCVEMPulses','OfflineIceTopSLCVEMPulses'],
  Output='IceTopVEMPulses',
  streams=[icetray.I3Frame.DAQ],
  )
tray.AddModule(checkDeltaT, "dT",
				# streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
				)



# tray.AddModule("I3Writer","i3writer",
# 	filename=str(outputDir)+dataSetClean+str(fileName)+"CleanVEMEvts.i3.gz",
# 	streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
# 						)
tray.Execute()
tray.Finish()