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

import numpy as np

from icecube.weighting.fluxes import GaisserH4a_IT
from icecube.weighting.weighting import icetop_mc_weights

weight_file = "/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json"
# weight_file = "/home/enpaudel/icecube/triggerStudy/simFiles/dataset_info_part.json" #my sim json file
#bash: python triggerProc.py /data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/12360/Level3_IC86.2012_12360_Run000001.i3.gz
#python /data/sim/IceTop/2012/filtered/CORSIKA-ice-top/12360/level2/0000000-0000999/Level2_IC86_corsika_icetop.010410.000050.i3.bz2 


exceptionTanks_HG = {39:62,26:62,67:64}
exceptionTanks_LG = {26:61,67:63}


CORSIKA_ID = "DAT059871"
outputDir = "/home/enpaudel/icecube/triggerStudy/simFiles/"

# GCD="/data/user/kath/testdata/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"
GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"
charge_threshold = 10**-3 #threshold for charges of afterpulse in units of vem
time_threshold = 10.0**6.0 #in units of ns

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/ITGenTest/FeDAT000001GenDet.i3.bz2", help='Input files after running detector.py.')
args = parser.parse_args()

fileName = args.input[0].split("/")[-1]
# fileName = fileName.split(".")[1] #for level3
fileName = str(fileName.split(".")[1])+"_"+ str(fileName.split(".")[2]) #for level2
print(fileName)
print("inputs",args.input)


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
		if frame.Has("QTriggerHierarchy"):
			triggerHierarchy = frame["QTriggerHierarchy"]
		elif frame.Has("I3TriggerHierarchy"):
			triggerHierarchy = frame["I3TriggerHierarchy"]
		# print(frame.keys())
		self.totalEvents += 1
		print("total no of events",self.totalEvents)
		if len(triggerHierarchy) == 0:
			frame["ITSMTTriggered"] = dataclasses.I3Double(0)
			frame["ITGlobalTriggered"] = dataclasses.I3Double(0)
			self.GlobalUntriggered+=1
			self.SMTUntriggered+=1
		else:
			SMTTriggers = [t for t in triggerHierarchy if (t.key.config_id == 102 and t.fired)]
			GlobalTriggers = [t for t in triggerHierarchy if t.key.config_id == None and t.fired]
			if len(SMTTriggers) != 0:
				frame["ITSMTTriggered"] = dataclasses.I3Double(1)
				self.SMTTriggered+=1
			else:
				frame["ITSMTTriggered"] = dataclasses.I3Double(0)
				self.SMTUntriggered+=1
			if len(GlobalTriggers) != 0:
				frame["ITGlobalTriggered"] = dataclasses.I3Double(1)
				self.GlobalTriggered+=1
			else:
				frame["ITGlobalTriggered"] = dataclasses.I3Double(0)
				self.GlobalUntriggered+=1
		self.PushFrame(frame)

	def Finish(self):
		print("total number of events",self.totalEvents)
		# print("number of SMT triggered events",self.SMTTriggered )
		# print("number of SMT untriggered events",self.SMTUntriggered )
		# print("number of Global triggered events",self.GlobalTriggered)
		# print("number of Gloabl untriggered events",self.GlobalUntriggered)
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

def timeCleanedPulses(frame,keys):
	for ikey in keys:
		if frame.Has(ikey):
			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,ikey)
			cleanEvt = []
			psmClean = dataclasses.I3RecoPulseSeriesMap()
			for omkey,pulses in psm:
				ps = dataclasses.I3RecoPulseSeries()
				for pulse in pulses:
					if pulse.time < time_threshold:
						ps.append(pulse)
				if len(ps) > 0:
					psmClean[omkey] = dataclasses.I3RecoPulseSeries(sorted(ps,
	            key=lambda pulse: pulse.time))
			frame[ikey+"CleanTime"] = psmClean

def chargeCleanedPulses(frame,keys):
	"""remove pulses that are below 0.16 threshold"""
	for ikey in keys:
		if frame.Has(ikey):
			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,ikey)
			cleanEvt = []
			psmClean = dataclasses.I3RecoPulseSeriesMap()
			for omkey,pulses in psm:
				ps = dataclasses.I3RecoPulseSeries()
				for pulse in pulses:
					if pulse.charge > charge_threshold:
						ps.append(pulse)
				if len(ps) > 0:
					psmClean[omkey] = dataclasses.I3RecoPulseSeries(sorted(ps,
	            key=lambda pulse: pulse.time))
			frame[ikey+"CleanCharge"] = psmClean

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
		# if NCh == 0:
		# 	# print("psm",psm)
		frame[str(pulseseries)+"TotalHit"] = dataclasses.I3Double(NCh)
		# print("No. of Hits in ",pulseseries,NCh,len(channels),frame[str(pulseseries)+"TotalTankHit"])

def deltaTHLCHit(frame,pulseseriesList):
	'''calculates duration of HLC hits;
	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	'''
	# if int(frame["I3EventHeader"].event_id) == int(8351):
	for pulseseries in pulseseriesList:
		psm = frame[pulseseries]
		if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
		times_list = []
		hit_stations = []
		hit_omkeys = []
		for om,pulses in psm:
			for pulse in pulses:
				hit_stations.append(om[0])
				hit_omkeys.append(om)
				break
		hit_stations = list(set(hit_stations))
		hit_tanks = []
		for istation in hit_stations:
			tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
			(iom[0] == istation and iom[0] in [26,39] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
			pulses = [psm[om] for om in tanks]
			hit_tanks += tanks
			# ipulses = [ipulse for ipulse in pulse for pulse in pulses]
			hit_times = [pulse[0].time for pulse in pulses]
			# hit_times = [ipulse.time for pulse in pulses for ipulse in pulse]
			if len(hit_times)>=2:
				deltaT = abs(hit_times[1]-hit_times[0])
				times_list.append(deltaT)
		frame[str(pulseseries)+"_delta_t"] = dataclasses.I3VectorFloat(times_list)
		frame[str(pulseseries)+"_hitTanks"] = dataclasses.I3Double(len(list(set(hit_tanks))))
		frame[str(pulseseries)+"_hitStations"] = dataclasses.I3Double(len(hit_stations))
		if len(frame[str(pulseseries)+"_delta_t"])>0 and min(frame[str(pulseseries)+"_delta_t"]) < 1000:
			frame[str(pulseseries)+"_isSTA1"] = dataclasses.I3Double(1)
		else:
			frame[str(pulseseries)+"_isSTA1"] = dataclasses.I3Double(0)

def checkFilter(frame,filterKeys):
	for filterKey in filterKeys:
		if frame.Has("QFilterMask"):
			filters = frame["QFilterMask"]
			print("filters",filters)
			print("writing in frame")
			print(filters[filterKey].condition_passed and filters[filterKey].prescale_passed)
			print(filters[filterKey].condition_passed and filters[filterKey].prescale_passed)
			print(filterKey+"_filter")
			if filters[filterKey].condition_passed and filters[filterKey].prescale_passed:
				frame[filterKey+"_filter"] = dataclasses.I3Double(1)
			else:
				frame[filterKey+"_filter"] = dataclasses.I3Double(0)
			print(frame[filterKey+"_filter"])



def deltaTSLCHit(frame,pulseseriesList):
	'''calculates duration of HLC hits;
	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	'''
	# if int(frame["I3EventHeader"].event_id) == int(8351):
	for pulseseries in pulseseriesList:
		psm = frame[pulseseries]
		if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
		times_list = []
		hit_stations = []
		hit_omkeys = []
		for om,pulses in psm:
			for pulse in pulses:
				hit_stations.append(om[0])
				hit_omkeys.append(om)
				break
		hit_stations = list(set(hit_stations))
		hit_tanks = []
		for istation in hit_stations:
			tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
			(iom[0] == istation and iom[0] in [26,39] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
			pulses = [psm[om] for om in tanks]
			hit_tanks += tanks
			# ipulses = [ipulse for ipulse in pulse for pulse in pulses]
			for pulse in pulses:
				times_list.append(pulse[0].time)
		times_list.sort()
		frame[str(pulseseries)+"_delta_t"] = dataclasses.I3VectorFloat(np.diff(times_list))
		frame[str(pulseseries)+"_hitTanks"] = dataclasses.I3Double(len(list(set(hit_tanks))))
		frame[str(pulseseries)+"_hitStations"] = dataclasses.I3Double(len(hit_stations))


def AddTimeSLCHLCTankHit(frame,pulseseriesList):
	'''calculates total SLC or HLC duration
	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	'''
	# print("Adding Hit duration")
	for pulseseries in pulseseriesList:
		timeList = []
		psm = frame[pulseseries]
		if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
		for om,pulses in psm:
			for pulse in pulses:
				timeList.append(pulse.time)
				# print("time check",om,pulse.time)
				break
		# print("timeList",pulseseries,timeList)
		if len(timeList) == 0:
			frame[str(pulseseries)+"HitTimeDuration"] = dataclasses.I3Double(-1)
		elif len(timeList) == 1:
			frame[str(pulseseries)+"HitTimeDuration"] = dataclasses.I3Double(-1*timeList[0])
		else:
			frame[str(pulseseries)+"HitTimeDuration"] = dataclasses.I3Double(max(timeList)-min(timeList))
			if frame[str(pulseseries)+"HitTimeDuration"] > time_threshold:
				frame[str(pulseseries)+"UnusualTime"] = dataclasses.I3Double(max(timeList)-min(timeList))
				print("Yell that there is  a problem", frame[str(pulseseries)+"HitTimeDuration"],max(timeList),min(timeList),"in file",fileName,"event",frame["I3EventHeader"],frame["I3EventHeader"].event_id)
		# print("time Duration",frame[str(pulseseries)+"HitTimeDuration"])
		# print("No. of Hits in ",pulseseries,NCh,len(channels),frame[str(pulseseries)+"TotalTankHit"])


def printKeyInFrames(frame):
	for ikey in frame.keys():
		print(ikey)

# def calcWeight(frame):
# 	flux = GaisserH4a_IT()
# 	p_energy = frame["MCPrimary"].energy 
# 	p_type = frame["MCPrimary"].type
# 	p_zenith = frame["MCPrimary"].dir.zenith
# 	print("trying to debug",p_energy,type(p_energy),p_type,p_zenith)
# 	if str(p_type) == "PPlus" and float(p_energy) < 10.0**7.0:
# 		dset = 123600
# 	elif str(p_type) == "PPlus" and float(p_energy) >= 10.0**7.0:
# 		dset = 123601
# 	elif str(p_type) == "He4Nucleus" and float(p_energy) < 10.0**7.0:
# 		dset = 126300
# 	elif str(p_type) == "He4Nucleus" and float(p_energy) >= 10.0**7.0:
# 		dset = 126301
# 	elif str(p_type) == "O16Nucleus" and float(p_energy) < 10.0**7.0:
# 		dset = 126310
# 	elif str(p_type) == "O16Nucleus" and float(p_energy) >= 10.0**7.0:
# 		dset = 126311
# 	elif str(p_type) == "Fe56Nucleus" and float(p_energy) < 10.0**7.0 :
# 		dset = 123620
# 	elif str(p_type) == "Fe56Nucleus" and float(p_energy) >= 10.0**7.0 :
# 		dset = 123621
# 	print("dset",dset,weight_file)
# 	generator = icetop_mc_weights(dset,dataset_file=weight_file)
# 	p_weights = flux(p_energy,p_type)/generator(p_energy,p_type,np.cos(p_zenith))
# 	frame["H4aWeight"]=dataclasses.I3Double(p_weights)

def calcWeight(frame):
	flux = GaisserH4a_IT()
	p_energy = frame["MCPrimary"].energy 
	p_type = frame["MCPrimary"].type
	p_zenith = frame["MCPrimary"].dir.zenith
	print("trying to debug",p_energy,type(p_energy),p_type,p_zenith)
	if str(p_type) == "PPlus":
		dset = 12360
	elif str(p_type) == "He4Nucleus":
		dset = 12630
	elif str(p_type) == "O16Nucleus":
		dset = 12631
	elif str(p_type) == "Fe56Nucleus":
		dset = 12362
	print("dset",dset,weight_file)
	generator = icetop_mc_weights(dset,dataset_file=weight_file)
	p_weights = flux(p_energy,p_type)/generator(p_energy,p_type,np.cos(p_zenith))
	frame["H4aWeight"]=dataclasses.I3Double(p_weights)


tray = I3Tray()
tray.AddModule("I3Reader","reader",
	           filenameList=[GCD]+args.input,
	           # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )
# tray.AddModule(TriggerCheck, "trigChk",
# 				# Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
# 				)
tray.AddModule(checkFilter,"filterChk",	
				# filterKeys=["IceTopSTA5_13","SDST_IceTopSTA3_13"],
				filterKeys=["IceTopSTA5_12","SDST_IceTopSTA3_12"],
				Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
				)

tray.AddModule('I3HLCTankPulseMerger',
	           # filenameList=[GCD]+args.input,
	           InputVEMPulses = 'OfflineIceTopSLCVEMPulses',
	           OutputTankPulses = 'OfflineIceTopSLCTankPulses',
	           ExcludedTanks  = 'ExcludedSLCTanks'
	           )

# tray.AddModule(timeCleanedPulses,"cleanTime",
# 	keys=['OfflineIceTopSLCVEMPulses','OfflineIceTopHLCVEMPulses',"IceTopPulses"],
# 	streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

# tray.AddModule(chargeCleanedPulses,"cleanCharge",
# 	keys=['OfflineIceTopSLCVEMPulsesCleanTime','OfflineIceTopHLCVEMPulsesCleanTime',"IceTopPulsesCleanTime"],
# 	streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

# tray.AddModule('I3HLCTankPulseMerger',"slcMerger2",
# 	           InputVEMPulses ='OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge',
# 	           OutputTankPulses = 'OfflineIceTopSLCTankPulsesCleanTimeCleanCharge',
# 	           ExcludedTanks  = 'ExcludedSLCTanksAfterCleaning'
# 	           )
# tray.AddModule('I3HLCTankPulseMerger',"hlcMerger",
# 	           InputVEMPulses ='OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge',
# 	           OutputTankPulses = 'OfflineIceTopHLCTankPulsesCleanTimeCleanCharge',
# 	           ExcludedTanks  = 'ExcludedHLCTanksAfterCleaning'
# 	           )

tray.AddModule(AddTotalCharge,"addCharge",
	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           keys=['OfflineIceTopSLCTankPulses','OfflineIceTopHLCTankPulses',
	           		]
	           	)
tray.AddModule(AddTotalTankHit,"addHit",
	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           pulseseriesList=['OfflineIceTopSLCTankPulses','OfflineIceTopHLCTankPulses',
	           		]
	           # pulseseriesList=['OfflineIceTopHLCTankPulses']
	           	)
tray.AddModule(AddTimeSLCHLCTankHit,"addTime",
	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           pulseseriesList=['OfflineIceTopSLCTankPulses','OfflineIceTopHLCTankPulses',]
	           # pulseseriesList=['OfflineIceTopHLCTankPulses']
	           	)

tray.AddModule(deltaTHLCHit,"deltT",
	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           # Streams=[icetray.I3Frame.DAQ],
	           pulseseriesList=['OfflineIceTopHLCTankPulses'
	           ]
	           # pulseseriesList=['OfflineIceTopHLCTankPulses']
	           	)
tray.AddModule(deltaTSLCHit,"deltTSLC",
	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           # Streams=[icetray.I3Frame.DAQ],
	           pulseseriesList=['OfflineIceTopSLCTankPulses'
	           ]
	           # pulseseriesList=['OfflineIceTopHLCTankPulses']
	           	)
#########################################
#########################################


# keep_list = ["BeaconLaunches","I3Triggers","MCPrimary","ClusterCleaningExcludedTanks","OfflineIceTopHLCTankPulses",
# 			"FilterMask_NullSplit0","PassedKeepSuperDSTOnly","CalibratedWaveformRange","I3EventHeader",
# 			"OfflineIceTopHLCTankPulsesHitTimeDuration","InIceDSTPulses","PassedAnyFilter","I3MCPESeriesMap",
# 			"I3Calibration","IceTopComponentPulses_Electron","I3VEMCalData","DSTTriggers","I3MCPulseSeriesMapParticleIDMap",
# 			"BadDomsList","I3MCTree","I3SuperDST","QFilterMask","I3MCPulseSeriesMap","IceTopComponentPulses_ElectronFromChargedMesons",
# 			"OfflineIceTopSLCVEMPulses","IceTopComponentPulses_Gamma","IceTopComponentPulses_GammaFromChargedMesons",
# 			"TankPulseMergerExcludedTanks","IceTopComponentPulses_Muon","IceTopDSTPulses","IceTopPulses","IceTopRawData",
# 			"I3MCPESeriesMapParticleIDMap","InIceRawData","PassedConventional","OfflineIceTopSLCVEMPulsesTotalHit",
# 			"ExcludedSLCTanks","UncleanedInIcePulsesTimeRange","MCPrimaryInfo","CleanIceTopRawData","OfflineIceTopHLCVEMPulses",
# 			"OfflineIceTopHLCVEMPulsesTotalCharge","IceTopComponentPulses_Hadron","InIcePulses","SimTrimmer","QTriggerHierarchy"]
# keep_list = ["I3Triggers","MCPrimary","ClusterCleaningExcludedTanks","OfflineIceTopHLCTankPulses",
# 			"FilterMask_NullSplit0","PassedKeepSuperDSTOnly","CalibratedWaveformRange","I3EventHeader",
# 			"OfflineIceTopHLCTankPulsesHitTimeDuration","InIceDSTPulses","PassedAnyFilter","I3MCPESeriesMap",
# 			"I3Calibration","IceTopComponentPulses_Electron","I3VEMCalData","DSTTriggers","I3MCPulseSeriesMapParticleIDMap",
# 			"BadDomsList","I3MCTree","I3SuperDST","QFilterMask","I3MCPulseSeriesMap","IceTopComponentPulses_ElectronFromChargedMesons",
# 			"OfflineIceTopSLCVEMPulses","IceTopComponentPulses_Gamma","IceTopComponentPulses_GammaFromChargedMesons",
# 			"TankPulseMergerExcludedTanks","IceTopComponentPulses_Muon","IceTopDSTPulses","IceTopPulses","IceTopRawData",
# 			"I3MCPESeriesMapParticleIDMap","InIceRawData","PassedConventional","OfflineIceTopSLCVEMPulsesTotalHit",
# 			"ExcludedSLCTanks","UncleanedInIcePulsesTimeRange","MCPrimaryInfo","CleanIceTopRawData","OfflineIceTopHLCVEMPulses",
# 			"OfflineIceTopHLCVEMPulsesTotalCharge","IceTopComponentPulses_Hadron","InIcePulses","SimTrimmer","QTriggerHierarchy",
# 			"MCPrimary","ITSMTTriggered","ITGlobalTriggered","OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulsesTotalCharge",
# 			"OfflineIceTopHLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalHit","OfflineIceTopSLCTankPulsesTotalHit","OfflineIceTopHLCVEMPulses",
# 			"OfflineIceTopSLCVEMPulses","OfflineIceTopSLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalHit",
# 			"OfflineIceTopSLCVEMPulsesTotalHit","OfflineIceTopSLCTankPulsesHitTimeDuration","OfflineIceTopHLCTankPulsesHitTimeDuration","UnusualTime"
# 			]

keep_list = ["I3Triggers","MCPrimary","OfflineIceTopHLCTankPulses","I3EventHeader",
			"OfflineIceTopHLCTankPulsesHitTimeDuration","InIceDSTPulses","PassedAnyFilter","I3MCPESeriesMap",
			"I3Calibration","IceTopComponentPulses_Electron","I3VEMCalData","DSTTriggers","I3MCPulseSeriesMapParticleIDMap",
			"BadDomsList","I3MCTree","I3SuperDST","QFilterMask","I3MCPulseSeriesMap","IceTopComponentPulses_ElectronFromChargedMesons",
			"OfflineIceTopSLCVEMPulses","IceTopComponentPulses_Gamma","IceTopComponentPulses_GammaFromChargedMesons",
			"IceTopComponentPulses_Muon","IceTopDSTPulses","IceTopPulses","IceTopRawData","InIceRawData","OfflineIceTopSLCVEMPulsesTotalHit",
			"UncleanedInIcePulsesTimeRange","MCPrimaryInfo","CleanIceTopRawData","OfflineIceTopHLCVEMPulses",
			"OfflineIceTopHLCVEMPulsesTotalCharge","IceTopComponentPulses_Hadron","InIcePulses","SimTrimmer","QTriggerHierarchy",
			"MCPrimary","ITSMTTriggered","ITGlobalTriggered","OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulsesTotalCharge",
			"OfflineIceTopHLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalHit","OfflineIceTopSLCTankPulsesTotalHit","OfflineIceTopHLCVEMPulses",
			"OfflineIceTopSLCVEMPulses","OfflineIceTopSLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalHit",
			"OfflineIceTopSLCVEMPulsesTotalHit","OfflineIceTopSLCTankPulsesHitTimeDuration","OfflineIceTopHLCTankPulsesHitTimeDuration","UnusualTime"
			]
delList = ["BeaconLaunches","CalibratedWaveformRange","ClusterCleaningExcludedTanks","DSTTriggers","ExcludedSLCTanks","ExcludedSLCTanksAfterCleaning",
			"FilterMask_NullSplit0","I3DSTHeader","PassedKeepSuperDSTOnly","SimTrimmer","TankPulseMergerExcludedTanks","MCPrimaryInfo",
			"UncleanedInIcePulsesTimeRange","PassedConventional","PassedAnyFilter","I3MCPESeriesMapParticleIDMap"
			]
def delUnwantedFrameKeys(frame,keys):
	"""delete unwanted frame keys in list keys
	"""
	for ikey in keys:
		if frame.Has(ikey):
			del frame[ikey]

tray.AddModule(delUnwantedFrameKeys,'del_frames',
				keys=delList,
				streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
				)

# tray.AddModule(printKeyInFrames,'keys',
# 				streams=[icetray.I3Frame.DAQ]
# 				)
# def haveUnusualTime(frame):
# 	'''acts as  a simple filter
# 	'''
# 	if frame.Has("UnusualTime"):
# 		return True
# 	else:
# 		return False

# def haveHadronEvent(frame):
# 	'''acts as  a simple filter
# 	'''
# 	if frame.Has("IceTopComponentPulses_Hadron"):
# 		psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,"IceTopComponentPulses_Hadron")
# 		hadronEvt = []
# 		triggeredHadronEvt = []
# 		for omkey,pulses in psm:
# 			for pulse in pulses:
# 				if abs(pulse.charge - 0) > 0.1:
# 					hadronEvt.append(omkey)
# 					print("this event has hadron with key",omkey,pulse.charge,pulse.time)
# 					if frame["ITSMTTriggered"]==dataclasses.I3Double(0):
# 						print("trigger status",frame["ITSMTTriggered"])
# 						triggeredHadronEvt.append(omkey)
# 		if len(triggeredHadronEvt) > 0:
# 			return True
# 		else:
# 			return False
# 	else:
# 		return False


def haveHadronEvent(frame):
	'''acts as  a simple filter
	'''
	if frame.Has("IceTopComponentPulses_Hadron"):
		psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,"IceTopComponentPulses_Hadron")
		hadronEvt = []
		for omkey,pulses in psm:
			for pulse in pulses:
				if abs(pulse.charge - 0) > 0.1:
					hadronEvt.append(omkey)
					print("this event has hadron with key",omkey,pulse.charge,pulse.time)
		if len(hadronEvt) > 0:
			return True
		else:
			return False
	else:
		return False

def haveElectronEvent(frame):
	'''acts as  a simple filter
	'''
	if frame.Has("IceTopComponentPulses_Electron"):
		psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,"IceTopComponentPulses_Electron")
		hadronEvt = []
		for omkey,pulses in psm:
			for pulse in pulses:
				if abs(pulse.charge - 0) > 0.1:
					hadronEvt.append(omkey)
					print("this event has hadron with key",omkey,pulse.charge,pulse.time)
		if len(hadronEvt) > 0:
			return True
		else:
			return False
	else:
		return False

def haveMuonEvent(frame):
	'''acts as  a simple filter
	'''
	if frame.Has("IceTopComponentPulses_Muon"):
		psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,"IceTopComponentPulses_Muon")
		hadronEvt = []
		for omkey,pulses in psm:
			for pulse in pulses:
				if abs(pulse.charge - 0) > 0.1:
					hadronEvt.append(omkey)
					print("this event has hadron with key",omkey,pulse.charge,pulse.time)
		if len(hadronEvt) > 0:
			return True
		else:
			return False
	else:
		return False

def nonHadronEvent(frame):
	'''acts as  a simple filter
	'''
	if frame.Has("IceTopComponentPulses_Hadron"):
		psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,"IceTopComponentPulses_Hadron")
		hadronEvt = []
		for omkey,pulses in psm:
			for pulse in pulses:
				if abs(pulse.charge - 0) > 0.1:
					hadronEvt.append(omkey)
					print("this event has hadron with key",omkey,pulse.charge,pulse.time)
		if len(hadronEvt) > 0:
			return False
		else:
			return True
	else:
		return True



		
def onlyTriggeredEvent(frame):
	'''acts as  a simple filter
	'''
	if frame["ITSMTTriggered"]==dataclasses.I3Double(1):
		return True
	else:
		return False

def onlyNonTriggeredEvent(frame):
	'''acts as  a simple filter
	'''
	if frame["ITSMTTriggered"]==dataclasses.I3Double(0):
		return True
	else:
		return False
def onlyProperTimeEvents(frame):
	'''acts as  a simple filter
	'''
	if frame["OfflineIceTopSLCTankPulsesHitTimeDuration"] < time_threshold:
		return True
	else:
		return False

def onlyImProperTimeEvents(frame):
	'''acts as  a simple filter
	'''
	if frame["OfflineIceTopSLCTankPulsesHitTimeDuration"] >= time_threshold:
		return True
	else:
		return False

def onlyProperTimeEventsHLC(frame):
	'''acts as  a simple filter
	'''
	if frame["OfflineIceTopHLCTankPulsesHitTimeDuration"] < time_threshold:
		return True
	else:
		return False

def onlyImProperTimeEventsHLC(frame):
	'''acts as  a simple filter
	'''
	if frame["OfflineIceTopHLCTankPulsesHitTimeDuration"] >= time_threshold:
		return True
	else:
		return False


# tray.AddModule(lambda frame : frame.Has("UnusualTime"),streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

# tray.AddModule("Keep","keep unusual time",
#  	           keys = keep_list,
#  	           If = haveUnusualTime,
#  	           )
# tray.AddModule(haveUnusualTime,"unusualTime",
# 				streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
# 				)
# tray.AddModule(haveHadronEvent,"hadronEvent",
# 				streams=[icetray.I3Frame.DAQ]
# 				)
# tray.AddModule(nonHadronEvent,"nonHadronEvent",
# 				streams=[icetray.I3Frame.DAQ]
# 				)
# tray.AddModule(onlyNonTriggeredEvent,"nonTrigEvent",
# 				streams=[icetray.I3Frame.DAQ]
# 				)
# tray.AddModule(onlyImProperTimeEvents,"ImpropTime",
# 				streams=[icetray.I3Frame.DAQ]
				# )
# tray.AddModule(onlyTriggeredEvent,"TrigEvent",
# 				streams=[icetray.I3Frame.DAQ]
# 				)
# tray.AddModule(onlyImProperTimeEventsHLC,"ImpropTimeHLC",
# 				streams=[icetray.I3Frame.DAQ]
# 				)
# tray.AddModule(onlyProperTimeEvents,"propTime",
# 				streams=[icetray.I3Frame.DAQ]
 				# )

tray.AddModule(calcWeight,"calcWeight",
				streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
 				)


tray.AddModule("I3Writer","i3writer",
	          # filename=str(outputDir)+"/hadronTimeTest/"+str(fileName)+"HadronNonTrigImProperTEvts.i3.gz",
	          # filename=str(outputDir)+"/dataSetCleanTest/"+str(fileName)+"CleanVEMEvts.i3.gz",
	          filename=str(outputDir)+"/dataSetCleanOfficialL2Test/"+str(fileName)+"CleanVEMEvts.i3.gz",
	          streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	          # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	          # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )

# tray.AddModule("I3NullSplitter","nullsplitter")

tray.Add(hdfwriter.I3HDFWriter, 'hdfNull',
    # Output=str(outputDir)+"/dataSetCleanTest/"+str(fileName)+"CleanVEMEvts.hdf5",
    Output=str(outputDir)+"/dataSetCleanOfficialL2Test/"+str(fileName)+"CleanVEMEvts.hdf5",
    # Output=str(outputDir)+"/hadronTimeTest/"+str(fileName)+"NullHadronNonTrigImProperTEvts.hdf5",
    CompressionLevel=9,
    # SubEventStreams=['IceTopSplit'],
    # SubEventStreams=['NullSplit'],
    # SubEventStreams=["nullsplitter",'IceTopSplit',"nullsplitter",'NullSplit',]
    SubEventStreams=["ice_top"],
    # Streams=[icetray.I3Frame.DAQ],
    keys = [
    "MCPrimary","I3EventHeader","ITSMTTriggered","ITGlobalTriggered","OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses",
    "OfflineIceTopSLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalHit",
    "OfflineIceTopSLCTankPulsesTotalHit","OfflineIceTopHLCVEMPulses","OfflineIceTopSLCVEMPulses","OfflineIceTopSLCVEMPulsesTotalCharge",
    "OfflineIceTopHLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalHit","OfflineIceTopSLCVEMPulsesTotalHit",
    "OfflineIceTopSLCTankPulsesHitTimeDuration","OfflineIceTopHLCTankPulsesHitTimeDuration","IceTopPulses",
    "IceTopPulsesCleanTime","IceTopPulsesCleanTimeCleanCharge",'OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge','OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge',
    'OfflineIceTopSLCTankPulsesCleanTimeCleanCharge','OfflineIceTopHLCTankPulsesCleanTimeCleanCharge',"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalCharge",
    "OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalCharge","OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalCharge",
    "OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalHit",
    "OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalHit",'OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration',
    'OfflineIceTopHLCTankPulsesCleanTimeCleanChargeHitTimeDuration',"OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalHit",
    "OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalHit",'OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeHitTimeDuration',
    'OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeHitTimeDuration',"OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t",
    "OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_delta_t","OfflineIceTopHLCVEMPulses_delta_t",
    "OfflineIceTopSLCVEMPulses_delta_t","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_hitTanks",
    "OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_hitTanks","OfflineIceTopHLCVEMPulses_hitTanks",
    "OfflineIceTopSLCVEMPulses_hitTanks","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_hitStations",
    "OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_hitStations","OfflineIceTopHLCVEMPulses_hitStations",
    "OfflineIceTopSLCVEMPulses_hitStations",'OfflineIceTopHLCVEMPulses_isSTA1','OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_isSTA1',"H4aWeight",
    "OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_delta_t",
    "OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_delta_t","OfflineIceTopHLCTankPulses_delta_t",
    "OfflineIceTopSLCTankPulses_delta_t","OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_hitTanks",
    "OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_hitTanks","OfflineIceTopHLCTankPulses_hitTanks",
    "OfflineIceTopSLCTankPulses_hitTanks","OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_hitStations",
    "OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_hitStations","OfflineIceTopHLCTankPulses_hitStations",
    "OfflineIceTopSLCTankPulses_hitStations",'OfflineIceTopHLCTankPulses_isSTA1','OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_isSTA1',
    "IceTopSTA5_13_filter","SDST_IceTopSTA3_13_filter","IceTopSTA5_12_filter","SDST_IceTopSTA3_12_filter"
    ]
    )

# tray.Add(hdfwriter.I3HDFWriter, 'hdfIce',
#     Output=str(outputDir)+"/hadronTimeTest/"+str(fileName)+"IceTopNonHadronNonTrigImProperTEvts.hdf5",
#     # Output=str(outputDir)+"/hadronTimeTest/"+str(fileName)+"IceTopHadronNonTrigImProperTEvts.hdf5",
#     CompressionLevel=9,
#     SubEventStreams=['IceTopSplit'],
#     # SubEventStreams=['NullSplit'],
#     # SubEventStreams=["nullsplitter",'IceTopSplit',"nullsplitter",'NullSplit',]
#     # SubEventStreams=["ice_top"],
#     # Streams=[icetray.I3Frame.DAQ],
#     keys = [
#     "MCPrimary","ITSMTTriggered","ITGlobalTriggered","OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulsesTotalCharge",
#     "OfflineIceTopHLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalHit","OfflineIceTopSLCTankPulsesTotalHit","OfflineIceTopHLCVEMPulses",
#     "OfflineIceTopSLCVEMPulses","OfflineIceTopSLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalHit",
#     "OfflineIceTopSLCVEMPulsesTotalHit","OfflineIceTopSLCTankPulsesHitTimeDuration","OfflineIceTopHLCTankPulsesHitTimeDuration"
#     ]
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