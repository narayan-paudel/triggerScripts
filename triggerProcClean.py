#!/usr/bin/env python3


import os
import sys
from I3Tray import I3Tray, I3Units
from icecube.simprod import segments
from icecube.simprod.util import simprodtray, arguments
from icecube.simprod.util import ReadI3Summary, WriteI3Summary
from icecube import icetray, dataclasses, dataio, icetop_Level3_scripts
# from icecube import radcube
from icecube.simprod.util.simprodtray import RunI3Tray
import argparse
import icecube.icetray
import math
from icecube import phys_services, sim_services
from icecube import tableio, hdfwriter
from icecube.simprod.util import PrintContext
from icecube import topeventcleaning, tpx
# from icecube import top_background_simulator


import numpy as np

from icecube.weighting.fluxes import GaisserH4a_IT
from icecube.weighting.weighting import icetop_mc_weights

# weight_file = "/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json"
# weight_file = "/home/enpaudel/icecube/triggerStudy/simFiles/dataset_info_partUneven.json" #my sim json file
# weight_file = "/home/enpaudel/icecube/triggerStudy/simFiles/dataset_info_100.json" #my sim json file
# weight_file = "/home/enpaudel/icecube/triggerStudy/simFiles/dataset_info_10.json" #my sim json file
weight_file = "/home/enpaudel/icecube/triggerStudy/simFiles/dataset_info200Uneven.json" #my sim json file

exceptionTanks_HG = {39:62,26:62,67:64}
exceptionTanks_LG = {26:61,67:63}

ConfigIDMap = {102:(6,5000,"HLC6"),
103:(6,5000,"tank6"),104:(7,5000,"tank7"),105:(8,5000,"tank8"),106:(9,5000,"tank9"),107:(10,5000,"tank10"),
113:(6,4000,"tank6"),114:(7,4000,"tank7"),115:(8,4000,"tank8"),116:(9,4000,"tank9"),117:(10,4000,"tank10"),
123:(6,3000,"tank6"),124:(7,3000,"tank7"),125:(8,3000,"tank8"),126:(9,3000,"tank9"),127:(10,3000,"tank10"),
133:(6,2000,"tank6"),134:(7,2000,"tank7"),135:(8,2000,"tank8"),136:(9,2000,"tank9"),137:(10,2000,"tank10")}
#config id mapped with threshold,timewindow and label
print("keys",ConfigIDMap.keys())


CORSIKA_ID = "DAT059871"
# outputDir = "/home/enpaudel/icecube/triggerStudy/simFiles/"
outputDir = "/home/enpaudel/icecube/triggerStudy/simFiles/"
dataSetClean = "/dataSetClean/"
# dataSetClean = "/dataSetCleanWFRT/"
# dataSetClean = "/dataSetCleanFRT/"

# GCD="/data/user/kath/testdata/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"
# GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"
# GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"
GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"

f = dataio.I3File(GCD)
geoFrame = f.pop_frame()
while f.more() and not "I3Geometry" in geoFrame:
		geoFrame = f.pop_frame()

charge_threshold = 10**-3 #threshold for charges of afterpulse in units of vem
time_threshold = 10.0**6.0 #in units of ns


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/ITGenTest/FeDAT000001GenDet.i3.bz2", help='Input files after running detector.py.')
args = parser.parse_args()

#python triggerProcClean.py ../simFiles/dataSetUnique/FeDAT005775GenDetFiltProcUnique.i3.gz


fileName = args.input[0].split("/")[-1]
fileName = fileName.split(".")[0]
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
		for ikey in ConfigIDMap.keys():
			frame = self.checkTriggers(frame,triggerHierarchy, ikey)
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

# def deltaTHLCHit(frame,pulseseriesList):
# 	'''calculates duration of HLC hits;
# 	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
# 	'''
# 	# if int(frame["I3EventHeader"].event_id) == int(8351):
# 	for pulseseries in pulseseriesList:
# 		psm = frame[pulseseries]
# 		if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
# 			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
# 		times_list = []
# 		hit_stations = []
# 		hit_omkeys = []
# 		for om,pulses in psm:
# 			for pulse in pulses:
# 				hit_stations.append(om[0])
# 				hit_omkeys.append(om)
# 				break
# 		hit_stations = list(set(hit_stations))
# 		hit_tanks = []
# 		for istation in hit_stations:
# 			tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
# 			(iom[0] == istation and iom[0] in [26,39] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
# 			pulses = [psm[om] for om in tanks]
# 			hit_tanks += tanks
# 			# ipulses = [ipulse for ipulse in pulse for pulse in pulses]
# 			hit_times = [pulse[0].time for pulse in pulses]
# 			# hit_times = [ipulse.time for pulse in pulses for ipulse in pulse]
# 			print("length of hit times",len(hit_times))
# 			if len(hit_times)>=2:
# 				deltaT = abs(hit_times[1]-hit_times[0])
# 				times_list.append(deltaT)
# 		frame[str(pulseseries)+"_delta_t"] = dataclasses.I3VectorFloat(times_list)
# 		frame[str(pulseseries)+"_hitTanks"] = dataclasses.I3Double(len(list(set(hit_tanks))))
# 		frame[str(pulseseries)+"_hitStations"] = dataclasses.I3Double(len(hit_stations))
# 		if len(frame[str(pulseseries)+"_delta_t"])>0 and min(frame[str(pulseseries)+"_delta_t"]) < 1000:
# 			frame[str(pulseseries)+"_isSTA1"] = dataclasses.I3Double(1)
# 		else:
# 			frame[str(pulseseries)+"_isSTA1"] = dataclasses.I3Double(0)

# def getSPTime(particle,omkey,pulseTime):
# 	'''
# 	returns a shower plane time. Both script works exactly same
# 	'''
# 	domgeo = geoFrame["I3Geometry"].omgeo[omkey]
# 	x_om = domgeo.position.x
# 	y_om = domgeo.position.y
# 	z_om = domgeo.position.z
# 	t_offset = (particle.dir.x * (x_om - particle.pos.x) + \
#                                   particle.dir.y * (y_om - particle.pos.y) + \
#                                   particle.dir.z * (z_om - particle.pos.z))/dataclasses.I3Constants.c
# 	print("t_offset",t_offset+pulseTime)
# 	return pulseTime + t_offset



def getSPTime(particle,omkey,pulseTime):
	'''
	returns a shower plane time.
	'''
	domgeo = geoFrame["I3Geometry"].omgeo[omkey]
	# pos_sc = radcube.GetShowerFromIC(domgeo.position - particle.pos, particle.dir)
	# t_offset = pos_sc.z/dataclasses.I3Constants.c
	# print("deltaT",pulseTime + t_offset)
	# return pulseTime + t_offset
	return pulseTime

def deltaTHLCHit(frame,pulseseriesList):
	'''calculates duration of HLC hits and also puts trigger condition it will pass;
	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	1000 ns as per definition of HLC triggers
	'''
	# if int(frame["I3EventHeader"].event_id) == int(8351):
	particle = frame["MCPrimary"]
	for pulseseries in pulseseriesList:
		psm = frame[pulseseries]
		if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
		times_list = []
		times_listSC = []
		hit_stations = []
		hit_omkeys = []
		domList = []
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
			doms = [om for om in tanks]
			hit_tanks += tanks
			domList += doms
			# ipulses = [ipulse for ipulse in pulse for pulse in pulses]
			hit_times = [pulse[0].time for pulse in pulses]
			hit_timesSC = [getSPTime(particle,iom,ipulse[0].time) for iom,ipulse in zip(doms,pulses)]
			# hit_times = [ipulse.time for pulse in pulses for ipulse in pulse]
			# print("length of hit times",len(hit_times))
			if len(hit_times)>=2:
				deltaT = abs(hit_times[1]-hit_times[0])
				times_list.append(deltaT)
			if len(hit_timesSC)>=2:
				deltaTSC = abs(hit_timesSC[1]-hit_timesSC[0])
				times_listSC.append(deltaTSC)
		if len(times_list) > 0:
			frame[str(pulseseries)+"_hasHLC"] = dataclasses.I3Double(1)
		else:
			frame[str(pulseseries)+"_hasHLC"] = dataclasses.I3Double(0)
		frame[str(pulseseries)+"_delta_t"] = dataclasses.I3VectorFloat(times_list)
		frame[str(pulseseries)+"_hitTanks"] = dataclasses.I3Double(len(list(set(hit_tanks))))
		frame[str(pulseseries)+"_hitStations"] = dataclasses.I3Double(len(hit_stations))
		if len(frame[str(pulseseries)+"_delta_t"])>0 and min(frame[str(pulseseries)+"_delta_t"]) < 1000:
			frame[str(pulseseries)+"_isSTA1"] = dataclasses.I3Double(1)
		else:
			frame[str(pulseseries)+"_isSTA1"] = dataclasses.I3Double(0)
		if len(frame[str(pulseseries)+"_delta_t"])>1 and min(frame[str(pulseseries)+"_delta_t"]) < 1000:
			frame[str(pulseseries)+"_isSTA2"] = dataclasses.I3Double(1)
		else:
			frame[str(pulseseries)+"_isSTA2"] = dataclasses.I3Double(0)
		if len(frame[str(pulseseries)+"_delta_t"])>2 and min(frame[str(pulseseries)+"_delta_t"]) < 1000:
			frame[str(pulseseries)+"_isSTA3"] = dataclasses.I3Double(1)
		else:
			frame[str(pulseseries)+"_isSTA3"] = dataclasses.I3Double(0)
		frame[str(pulseseries)+"_delta_tSC"] = dataclasses.I3VectorFloat(times_listSC)
		if len(frame[str(pulseseries)+"_delta_tSC"])>0 and min(frame[str(pulseseries)+"_delta_tSC"]) < 1000:
			frame[str(pulseseries)+"_isSTA1SC"] = dataclasses.I3Double(1)
		else:
			frame[str(pulseseries)+"_isSTA1SC"] = dataclasses.I3Double(0)
		if len(frame[str(pulseseries)+"_delta_tSC"])>1 and min(frame[str(pulseseries)+"_delta_tSC"]) < 1000:
			frame[str(pulseseries)+"_isSTA2SC"] = dataclasses.I3Double(1)
		else:
			frame[str(pulseseries)+"_isSTA2SC"] = dataclasses.I3Double(0)



def checkFilter(frame,filterKeys):
	for filterKey in filterKeys:
		if frame.Has("QFilterMask"):
			filters = frame["QFilterMask"]
			# print("filters",filters)
			# print("writing in frame")
			# print(filters[filterKey].condition_passed and filters[filterKey].prescale_passed)
			# print(filters[filterKey].condition_passed and filters[filterKey].prescale_passed)
			# print(filterKey+"_filter")
			if filters[filterKey].condition_passed and filters[filterKey].prescale_passed:
				frame[filterKey+"_filter"] = dataclasses.I3Double(1)
			else:
				frame[filterKey+"_filter"] = dataclasses.I3Double(0)
			print(frame[filterKey+"_filter"])

# def deltaTSLCHit(frame,pulseseriesList):
# 	'''calculates duration of HLC hits;
# 	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
# 	'''
# 	# if int(frame["I3EventHeader"].event_id) == int(8351):
# 	for pulseseries in pulseseriesList:
# 		psm = frame[pulseseries]
# 		if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
# 			psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
# 		times_list = []
# 		hit_stations = []
# 		hit_omkeys = []
# 		for om,pulses in psm:
# 			for pulse in pulses:
# 				hit_stations.append(om[0])
# 				hit_omkeys.append(om)
# 				break
# 		hit_stations = list(set(hit_stations))
# 		hit_tanks = []
# 		for istation in hit_stations:
# 			tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
# 			(iom[0] == istation and iom[0] in [26,39] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
# 			pulses = [psm[om] for om in tanks]
# 			hit_tanks += tanks
# 			# ipulses = [ipulse for ipulse in pulse for pulse in pulses]
# 			for pulse in pulses:
# 				times_list.append(pulse[0].time)
# 		times_list.sort()
# 		frame[str(pulseseries)+"_delta_t"] = dataclasses.I3VectorFloat(np.diff(times_list))
# 		frame[str(pulseseries)+"_hitTanks"] = dataclasses.I3Double(len(list(set(hit_tanks))))
# 		frame[str(pulseseries)+"_hitStations"] = dataclasses.I3Double(len(hit_stations))

def nHitTrigger(particle,omList,pulseList,nHit,window,showerPlane):
	"""looks if trigger condition is met for nHit in given window in ns"""
	if showerPlane == True:
		timeList = [getSPTime(particle,iom,ipulse[0].time) for iom,ipulse in zip(omList,pulseList)]
	elif showerPlane == False:
		timeList = [ipulse[0].time for ipulse in pulseList]
	diffs = []
	triggerFlag = 0
	if len(timeList) < nHit:
		triggerFlag = 0
	else:
		for i in range(len(timeList)-(nHit-1)):
			idiff = timeList[i+(nHit-1)]-timeList[i]
			diffs.append(idiff)
			if idiff <= window:
				triggerFlag = 1
				break
		if min(diffs)<=window:
			triggerFlag = 1
		else:
			triggerFlag = 0
	return triggerFlag

def SLCTrigger(pulseList,nHit,window):
	"""looks if trigger condition is met for nHit in given window in ns a different approach"""
	timeList = [ipulse[0].time for ipulse in pulseList]
	timeList = np.asarray(timeList)
	triggerFlag = 0
	for itime in timeList:
		inTimeHits = [1 for i in timeList if itime <= i < itime+window]
		if len(inTimeHits) >= nHit:
			triggerFlag = 1
			break
	return triggerFlag

def hasHLCHit(frame,HLCpulseseriesList):
	HLCpulseSeries = HLCpulseseriesList[0]
	psm = frame[HLCpulseSeries]
	if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
		psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, HLCpulseSeries)
	times_listHLC = []
	hit_stationsHLC = []
	hit_omkeysHLC = []
	for om,pulses in psm:
		for pulse in pulses:
			hit_stationsHLC.append(om[0])
			hit_omkeysHLC.append(om)
			break
	hit_stationsHLC = list(set(hit_stationsHLC))
	hit_tanksHLC = []
	for istation in hit_stationsHLC:
		tanks = [iom for iom in hit_omkeysHLC if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
		(iom[0] == istation and iom[0] in [26,39] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
		pulses = [psm[om] for om in tanks]
		hit_times = [pulse[0].time for pulse in pulses]
		if len(hit_times)>=2:
			deltaT = abs(hit_times[1]-hit_times[0])
			times_listHLC.append(deltaT)
	if len(times_listHLC) > 0:
		hasHLC = True
	else:
		hasHLC = False

def deltaT3SLCHit(frame,SLCpulseseriesList,HLCpulseseriesList):
	'''calculates duration of first and third SLC hits;
	keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	'''
	# if int(frame["I3EventHeader"].event_id) == int(8351):
	# for pulseseries in pulseseriesList:
	particle = frame["MCPrimary"]

	#now looking at SLC pulses
	pulseseries = SLCpulseseriesList[0]
	psm = frame[pulseseries]
	if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
		psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
	times_list = []
	hit_stations = []
	hit_omkeys = []
	# print("pulse series",frame["I3EventHeader"].event_id,SLCpulseseriesList[0])
	for om,pulses in psm:
		for pulse in pulses:
			print("OMKeys",om[0],om[1],om[2])
			hit_stations.append(om[0])
			hit_omkeys.append(om)
			break
	hit_stations = list(set(hit_stations))
	hit_tanks = []
	SLCpulsesList = []
	domList = []
	for istation in hit_stations:
		tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
		(iom[0] == istation and iom[0] in [26,39] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
		pulses = [psm[om] for om in tanks]
		doms = [om for om in tanks]
		SLCpulsesList += pulses
		domList += doms
		hit_tanks += tanks
		# ipulses = [ipulse for ipulse in pulse for pulse in pulses]
	domList = [om for om,_ in sorted(zip(domList, SLCpulsesList), key=lambda pair: pair[1][0].time)]
	SLCpulsesList.sort(key=lambda x:x[0].time)
	times_list = [ipulse[0].time for ipulse in SLCpulsesList]
	frame[str(pulseseries)+"_hitTanks"] = dataclasses.I3Double(len(list(set(hit_tanks))))
	frame[str(pulseseries)+"_hitStations"] = dataclasses.I3Double(len(hit_stations))
	if len(SLCpulsesList)>4:
		frame[str(pulseseries)+"_2SLCduration_t"] = dataclasses.I3VectorFloat(np.diff(times_list))
		frame[str(pulseseries)+"_3SLCduration_t"] = dataclasses.I3Double(SLCpulsesList[2][0].time - SLCpulsesList[0][0].time)
		frame[str(pulseseries)+"_4SLCduration_t"] = dataclasses.I3Double(SLCpulsesList[3][0].time - SLCpulsesList[0][0].time)
		frame[str(pulseseries)+"_5SLCduration_t"] = dataclasses.I3Double(SLCpulsesList[4][0].time - SLCpulsesList[0][0].time)
		frame[str(pulseseries)+"_isSLC3"] = dataclasses.I3Double(SLCTrigger(particle,domList,SLCpulsesList,nHit=3,window=6000))
		frame[str(pulseseries)+"_isSLC4"] = dataclasses.I3Double(SLCTrigger(particle,domList,SLCpulsesList,nHit=4,window=6000))
		frame[str(pulseseries)+"_isSLC5"] = dataclasses.I3Double(SLCTrigger(particle,domList,SLCpulsesList,nHit=5,window=6000))
	elif len(SLCpulsesList)>3:			
		frame[str(pulseseries)+"_2SLCduration_t"] = dataclasses.I3VectorFloat(np.diff(times_list))
		frame[str(pulseseries)+"_3SLCduration_t"] = dataclasses.I3Double(SLCpulsesList[2][0].time - SLCpulsesList[0][0].time)
		frame[str(pulseseries)+"_4SLCduration_t"] = dataclasses.I3Double(SLCpulsesList[3][0].time - SLCpulsesList[0][0].time)
		frame[str(pulseseries)+"_isSLC3"] = dataclasses.I3Double(SLCTrigger(particle,domList,SLCpulsesList,nHit=3,window=6000))
		frame[str(pulseseries)+"_isSLC4"] = dataclasses.I3Double(SLCTrigger(particle,domList,SLCpulsesList,nHit=4,window=6000))
	elif len(SLCpulsesList)>2:
		frame[str(pulseseries)+"_2SLCduration_t"] = dataclasses.I3VectorFloat(np.diff(times_list))
		frame[str(pulseseries)+"_3SLCduration_t"] = dataclasses.I3Double(SLCpulsesList[2][0].time - SLCpulsesList[0][0].time)
		frame[str(pulseseries)+"_isSLC3"] = dataclasses.I3Double(SLCTrigger(particle,domList,SLCpulsesList,nHit=3,window=6000))
	elif len(SLCpulsesList)>1:
		frame[str(pulseseries)+"_2SLCduration_t"] = dataclasses.I3VectorFloat(np.diff(times_list))
	for deltaTstr in [str(pulseseries)+"_2SLCduration_t",str(pulseseries)+"_3SLCduration_t",str(pulseseries)+"_4SLCduration_t",str(pulseseries)+"_5SLCduration_t"]:
		if not frame.Has(deltaTstr):
			frame[deltaTstr] = dataclasses.I3Double(-111)
	for trigStr in [str(pulseseries)+"_isSLC5",str(pulseseries)+"_isSLC4",str(pulseseries)+"_isSLC3"]:
		if not frame.Has(trigStr):
			frame[trigStr] = dataclasses.I3Double(0)

def getPulses(frame,pulseseries):
	'''
	Get pulses from pulse series
	'''
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
	pulsesList = []
	domList = []
	for istation in hit_stations:
		tanks = [iom for iom in hit_omkeys if (iom[0] == istation and iom[0] not in exceptionTanks_HG.keys()) and (iom[1] == 61 or iom[1]== 63)  or 
		(iom[0] == istation and iom[0] in [26,39] and (iom[1] == 62 or iom[1]== 63) ) or (iom[0] == istation and iom[0] == 67 and (iom[1] == 61 or iom[1] == 64) )]
		pulses = [psm[om] for om in tanks]
		doms = [om for om in tanks]
		pulsesList += pulses
		domList += doms
		hit_tanks += tanks
	return pulsesList,domList,hit_tanks,hit_stations


def triggerSLCHLCHits(frame,SLCpulseseriesList,HLCpulseseriesList,suffix):
	'''
	looks into SLC and HLC hits; keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	'''
	particle = frame["MCPrimary"]
	pulseseriesSLC = SLCpulseseriesList[0]
	pulseseriesHLC = HLCpulseseriesList[0]
	pulseListSLC,domListSLC,hit_tanksSLC,hit_stationsSLC = getPulses(frame,pulseseriesSLC)
	pulseListHLC,domListHLC,hit_tanksHLC,hit_stationsHLC = getPulses(frame,pulseseriesHLC)
	pulseList = [*pulseListHLC,*pulseListSLC]
	domList = [*domListHLC,*domListSLC]
	hit_tanks = [*hit_tanksHLC,*hit_tanksSLC]
	hit_stations = [*hit_stationsHLC,*hit_stationsSLC]
	hit_stations = list(set(hit_stations))
	domList = [om for om,_ in sorted(zip(domList, pulseList), key=lambda pair: pair[1][0].time)]
	pulseList.sort(key=lambda x:x[0].time)
	times_list = [ipulse[0].time for ipulse in pulseList]
	times_listSC = [getSPTime(particle,idom,ipulse[0].time) for ipulse,idom in zip(pulseList,domList)]
	frame["HLCSLC_hitTanks"+suffix+""] = dataclasses.I3Double(len(list(set(hit_tanks))))
	frame["HLCSLC_hitStations"+suffix+""] = dataclasses.I3Double(len(hit_stations))
	if len(pulseList) == 0:
		frame["HLCSLC"+suffix+"_isTank1"] = dataclasses.I3Double(0)
	elif len(pulseList) > 0:
		frame["HLCSLC"+suffix+"_isTank1"] = dataclasses.I3Double(1)
	if len(pulseList) > 1:
		frame["HLCSLC_2TankHit"+suffix+"_t"] = dataclasses.I3VectorFloat(np.diff(times_list))
		frame["HLCSLC_2TankHit"+suffix+"_tSC"] = dataclasses.I3VectorFloat(np.diff(times_listSC))
		frame["HLCSLC"+suffix+"_isTank2"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=2,window=6000,showerPlane=False))
		if len(pulseList) > 2:
			frame["HLCSLC_3TankHit"+suffix+"_t"] = dataclasses.I3Double(pulseList[2][0].time - pulseList[0][0].time)
			frame["HLCSLC_3TankHit"+suffix+"_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[2],pulseList[2][0].time)\
			 - getSPTime(particle,domList[0],pulseList[0][0].time))
			frame["HLCSLC"+suffix+"_isTank3"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=3,window=6000,showerPlane=False))
			if len(pulseList) > 3:
				frame["HLCSLC_4TankHit"+suffix+"_t"] = dataclasses.I3Double(pulseList[3][0].time - pulseList[0][0].time)
				frame["HLCSLC_4TankHit"+suffix+"_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[3],pulseList[3][0].time)\
				 - getSPTime(particle,domList[0],pulseList[0][0].time))
				frame["HLCSLC"+suffix+"_isTank4"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=4,window=6000,showerPlane=False))
				if len(pulseList) > 4:
					frame["HLCSLC_5TankHit"+suffix+"_t"] = dataclasses.I3Double(pulseList[4][0].time - pulseList[0][0].time)
					frame["HLCSLC_5TankHit"+suffix+"_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[4],pulseList[4][0].time)\
					 - getSPTime(particle,domList[0],pulseList[0][0].time))
					frame["HLCSLC"+suffix+"_isTank5"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=5,window=6000,showerPlane=False))
					if len(pulseList) > 5:
						frame["HLCSLC_6TankHit"+suffix+"_t"] = dataclasses.I3Double(pulseList[5][0].time - pulseList[0][0].time)
						frame["HLCSLC_6TankHit"+suffix+"_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[5],pulseList[5][0].time)\
						 - getSPTime(particle,domList[0],pulseList[0][0].time))
						frame["HLCSLC"+suffix+"_isTank6"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=6,window=6000,showerPlane=False))
						if len(pulseList) > 6:
							frame["HLCSLC_7TankHit"+suffix+"_t"] = dataclasses.I3Double(pulseList[6][0].time - pulseList[0][0].time)
							frame["HLCSLC_7TankHit"+suffix+"_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[6],pulseList[6][0].time)\
							 - getSPTime(particle,domList[0],pulseList[0][0].time))
							frame["HLCSLC"+suffix+"_isTank7"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=7,window=6000,showerPlane=False))
							if len(pulseList) > 7:
								frame["HLCSLC_8TankHit"+suffix+"_t"] = dataclasses.I3Double(pulseList[7][0].time - pulseList[0][0].time)
								frame["HLCSLC_8TankHit"+suffix+"_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[7],pulseList[7][0].time)\
								 - getSPTime(particle,domList[0],pulseList[0][0].time))
								frame["HLCSLC"+suffix+"_isTank8"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=8,window=6000,showerPlane=False))
	for deltaTstr in ["HLCSLC_2TankHit"+suffix+"_t","HLCSLC_2TankHit"+suffix+"_tSC"]:
		if not frame.Has(deltaTstr):
			frame[deltaTstr] = dataclasses.I3VectorFloat(np.diff(times_list))
	for deltaTstr in ["HLCSLC_3TankHit"+suffix+"_t","HLCSLC_4TankHit"+suffix+"_t","HLCSLC_5TankHit"+suffix+"_t","HLCSLC_6TankHit"+suffix+"_t","HLCSLC_7TankHit"+suffix+"_t","HLCSLC_8TankHit"+suffix+"_t","HLCSLC_3TankHit"+suffix+"_tSC","HLCSLC_4TankHit"+suffix+"_tSC",\
	"HLCSLC_5TankHit"+suffix+"_tSC","HLCSLC_6TankHit"+suffix+"_tSC","HLCSLC_7TankHit"+suffix+"_tSC","HLCSLC_8TankHit"+suffix+"_tSC"]:
		if not frame.Has(deltaTstr):
			frame[deltaTstr] = dataclasses.I3Double(-111)
	for trigStr in ["HLCSLC"+suffix+"_isTank5","HLCSLC"+suffix+"_isTank4","HLCSLC"+suffix+"_isTank3","HLCSLC"+suffix+"_isTank6","HLCSLC"+suffix+"_isTank7","HLCSLC"+suffix+"_isTank8"]:
		if not frame.Has(trigStr):
			frame[trigStr] = dataclasses.I3Double(0)

def triggerSLCHits(frame,SLCpulseseriesList):
	'''
	looks into SLC and HLC hits; keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
	'''
	particle = frame["MCPrimary"]
	pulseseries = SLCpulseseriesList[0]
	pulseList,domList,hit_tanks,hit_stations = getPulses(frame,pulseseries)
	domList = [om for om,_ in sorted(zip(domList, pulseList), key=lambda pair: pair[1][0].time)]
	pulseList.sort(key=lambda x:x[0].time)
	times_list = [ipulse[0].time for ipulse in pulseList]
	times_listSC = [getSPTime(particle,idom,ipulse[0].time) for ipulse,idom in zip(pulseList,domList)]
	frame["SLC_hitTanks"] = dataclasses.I3Double(len(list(set(hit_tanks))))
	frame["SLC_hitStations"] = dataclasses.I3Double(len(hit_stations))
	if len(pulseList) > 1:
		frame[str(pulseseries)+"_2SLCTankHit_t"] = dataclasses.I3VectorFloat(np.diff(times_list))
		frame[str(pulseseries)+"_2SLCTankHit_tSC"] = dataclasses.I3VectorFloat(np.diff(times_listSC))
		if len(pulseList) > 2:
			frame[str(pulseseries)+"_3SLCTankHit_t"] = dataclasses.I3Double(pulseList[2][0].time - pulseList[0][0].time)
			frame[str(pulseseries)+"_3SLCTankHit_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[2],pulseList[2][0].time)\
			 - getSPTime(particle,domList[0],pulseList[0][0].time))
			frame[str(pulseseries)+"_isSLCTank3"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=3,window=6000,showerPlane=False))
			if len(pulseList) > 3:
				frame[str(pulseseries)+"_4SLCTankHit_t"] = dataclasses.I3Double(pulseList[3][0].time - pulseList[0][0].time)
				frame[str(pulseseries)+"_4SLCTankHit_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[3],pulseList[3][0].time)\
				 - getSPTime(particle,domList[0],pulseList[0][0].time))
				frame[str(pulseseries)+"_isSLCTank4"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=4,window=6000,showerPlane=False))
				if len(pulseList) > 4:
					frame[str(pulseseries)+"_5SLCTankHit_t"] = dataclasses.I3Double(pulseList[4][0].time - pulseList[0][0].time)
					frame[str(pulseseries)+"_5SLCTankHit_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[4],pulseList[4][0].time)\
					 - getSPTime(particle,domList[0],pulseList[0][0].time))
					frame[str(pulseseries)+"_isSLCTank5"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=5,window=6000,showerPlane=False))
					if len(pulseList) > 5:
						frame[str(pulseseries)+"_6SLCTankHit_t"] = dataclasses.I3Double(pulseList[5][0].time - pulseList[0][0].time)
						frame[str(pulseseries)+"_6SLCTankHit_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[5],pulseList[5][0].time)\
						 - getSPTime(particle,domList[0],pulseList[0][0].time))
						frame[str(pulseseries)+"_isSLCTank6"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=6,window=6000,showerPlane=False))
						if len(pulseList) > 6:
							frame[str(pulseseries)+"_7SLCTankHit_t"] = dataclasses.I3Double(pulseList[6][0].time - pulseList[0][0].time)
							frame[str(pulseseries)+"_7SLCTankHit_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[6],pulseList[6][0].time)\
							 - getSPTime(particle,domList[0],pulseList[0][0].time))
							frame[str(pulseseries)+"_isSLCTank7"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=7,window=6000,showerPlane=False))
							if len(pulseList) > 7:
								frame[str(pulseseries)+"_8SLCTankHit_t"] = dataclasses.I3Double(pulseList[7][0].time - pulseList[0][0].time)
								frame[str(pulseseries)+"_8SLCTankHit_tSC"] = dataclasses.I3Double(getSPTime(particle,domList[7],pulseList[7][0].time)\
								 - getSPTime(particle,domList[0],pulseList[0][0].time))
								frame[str(pulseseries)+"_isSLCTank8"] = dataclasses.I3Double(nHitTrigger(particle,domList,pulseList,nHit=8,window=6000,showerPlane=False))

	for deltaTstr in [str(pulseseries)+"_2TankHit_t",str(pulseseries)+"_2TankHit_tSC"]:
		if not frame.Has(deltaTstr):
			frame[deltaTstr] = dataclasses.I3VectorFloat([])
	for deltaTstr in [str(pulseseries)+"_3SLCTankHit_t",str(pulseseries)+"_4SLCTankHit_t",str(pulseseries)+"_5SLCTankHit_t",
	str(pulseseries)+"_6SLCTankHit_t",str(pulseseries)+"_7SLCTankHit_t",str(pulseseries)+"_8SLCTankHit_t",
	str(pulseseries)+"_3SLCTankHit_tSC",str(pulseseries)+"_4SLCTankHit_tSC",str(pulseseries)+"_5SLCTankHit_tSC",
	str(pulseseries)+"_6SLCTankHit_tSC",str(pulseseries)+"_7SLCTankHit_tSC",str(pulseseries)+"_8SLCTankHit_tSC"]:
		if not frame.Has(deltaTstr):
			frame[deltaTstr] = dataclasses.I3Double(-111)
	for trigStr in [str(pulseseries)+"_isSLCTank5",str(pulseseries)+"_isSLCTank4",str(pulseseries)+"_isSLCTank3"\
	,str(pulseseries)+"_isSLCTank6",str(pulseseries)+"_isSLCTank7",str(pulseseries)+"_isSLCTank8"]:
		if not frame.Has(trigStr):
			frame[trigStr] = dataclasses.I3Double(0)




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
	# print("trying to debug",p_energy,type(p_energy),p_type,p_zenith)
	if str(p_type) == "PPlus":
		dset = 12360
	elif str(p_type) == "He4Nucleus":
		dset = 12630
	elif str(p_type) == "O16Nucleus":
		dset = 12631
	elif str(p_type) == "Fe56Nucleus":
		dset = 12362
	# print("dset",dset,weight_file)
	generator = icetop_mc_weights(dset,dataset_file=weight_file)
	p_weights = flux(p_energy,p_type)/generator(p_energy,p_type,np.cos(p_zenith))
	frame["H4aWeight"]=dataclasses.I3Double(p_weights)

name = ""


tray = I3Tray()
tray.AddModule("I3Reader","reader",
						 filenameList=[GCD]+args.input,
						 # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						 # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
						)
tray.AddModule(TriggerCheck, "trigChk",
				# Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
				)
tray.AddModule(checkFilter,"filterChk",	
				filterKeys=["IceTopSTA5_13","SDST_IceTopSTA3_13"],
				Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
				)

tray.AddModule('I3TankPulseMerger',
						 # filenameList=[GCD]+args.input,
						 InputVEMPulses = 'OfflineIceTopSLCVEMPulses',
						 OutputTankPulses = 'OfflineIceTopSLCTankPulses',
						 ExcludedTanks  = 'ExcludedSLCTanks'
						 )

tray.AddModule(timeCleanedPulses,"cleanTime",
	keys=['OfflineIceTopSLCVEMPulses','OfflineIceTopHLCVEMPulses',"IceTopPulses"],
	streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

tray.AddModule(chargeCleanedPulses,"cleanCharge",
	keys=['OfflineIceTopSLCVEMPulsesCleanTime','OfflineIceTopHLCVEMPulsesCleanTime',"IceTopPulsesCleanTime"],
	streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

tray.AddModule('I3TankPulseMerger',"slcMerger2",
						 InputVEMPulses ='OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge',
						 OutputTankPulses = 'OfflineIceTopSLCTankPulsesCleanTimeCleanCharge',
						 ExcludedTanks  = 'ExcludedSLCTanksAfterCleaning'
						 )
tray.AddModule('I3TankPulseMerger',"hlcMerger",
						 InputVEMPulses ='OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge',
						 OutputTankPulses = 'OfflineIceTopHLCTankPulsesCleanTimeCleanCharge',
						 ExcludedTanks  = 'ExcludedHLCTanksAfterCleaning'
						 )

# tray.AddModule('TopBackgroundSimulator', name + '_noisesim',
# 					 NoiseRate=1500.,#1424.
# 					 HLCTankPulses='OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge',
# 					 SLCTankPulses='OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge'
# 					 )

tray.AddModule(AddTotalCharge,"addCharge",
						 Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						 keys=['OfflineIceTopSLCTankPulses','OfflineIceTopHLCTankPulses','OfflineIceTopSLCVEMPulses','OfflineIceTopHLCVEMPulses',
								'OfflineIceTopSLCTankPulsesCleanTimeCleanCharge','OfflineIceTopHLCTankPulsesCleanTimeCleanCharge',
								'OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge','OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge']
							)
tray.AddModule(AddTotalTankHit,"addHit",
						 Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						 pulseseriesList=['OfflineIceTopSLCTankPulses','OfflineIceTopHLCTankPulses','OfflineIceTopSLCVEMPulses','OfflineIceTopHLCVEMPulses',
								'OfflineIceTopSLCTankPulsesCleanTimeCleanCharge','OfflineIceTopHLCTankPulsesCleanTimeCleanCharge',
								'OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge','OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge'
								]
						 # pulseseriesList=['OfflineIceTopHLCTankPulses']
							)
tray.AddModule(AddTimeSLCHLCTankHit,"addTime",
						 Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						 pulseseriesList=['OfflineIceTopSLCTankPulses','OfflineIceTopHLCTankPulses','OfflineIceTopSLCTankPulsesCleanTimeCleanCharge',
						 'OfflineIceTopHLCTankPulsesCleanTimeCleanCharge','OfflineIceTopSLCVEMPulses','OfflineIceTopHLCVEMPulses','OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge',
						 'OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge']
						 # pulseseriesList=['OfflineIceTopHLCTankPulses']
							)

tray.AddModule(deltaTHLCHit,"deltT",
						 Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						 # Streams=[icetray.I3Frame.DAQ],
						 pulseseriesList=['OfflineIceTopHLCVEMPulses','OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge',
						 'OfflineIceTopHLCTankPulses','OfflineIceTopHLCTankPulsesCleanTimeCleanCharge'
						 ]
						 # pulseseriesList=['OfflineIceTopHLCTankPulses']
							)
# tray.AddModule(deltaTSLCHit,"deltTSLC",
# 	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
# 	           # Streams=[icetray.I3Frame.DAQ],
# 	           pulseseriesList=['OfflineIceTopSLCVEMPulses','OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge',
# 	           'OfflineIceTopSLCTankPulses','OfflineIceTopSLCTankPulsesCleanTimeCleanCharge'
# 	           ]
# 	           # pulseseriesList=['OfflineIceTopHLCTankPulses']
# 	           	)

# tray.AddModule(deltaT3SLCHit,"deltT3SLC",
# 	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
# 	           # Streams=[icetray.I3Frame.DAQ],
# 	           SLCpulseseriesList=['OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge'
# 	           ],
# 	           HLCpulseseriesList=['OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge'
# 	           ]
# 	           # pulseseriesList=['OfflineIceTopHLCTankPulses']
# 	           	)


tray.AddModule(triggerSLCHLCHits,"deltTSLCHLC",
						 Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						 # Streams=[icetray.I3Frame.DAQ],
						 SLCpulseseriesList=['OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge'],
						 HLCpulseseriesList=['OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge'],
						 # pulseseriesList=['OfflineIceTopHLCTankPulses'],
						 suffix=""
							)
# tray.AddModule(triggerSLCHLCHits,"deltTSLCHLC",
# 	           Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
# 	           # Streams=[icetray.I3Frame.DAQ],
# 	           SLCpulseseriesList=['OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_NoBackground'],
# 	           HLCpulseseriesList=['OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_NoBackground'],
# 	           # pulseseriesList=['OfflineIceTopHLCTankPulses'],
# 	           suffix="NoBkg"
# 	           	)



tray.AddModule(triggerSLCHits,"deltTSLC",
						 Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						 # Streams=[icetray.I3Frame.DAQ],
						 SLCpulseseriesList=['OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge'
						 ],
						 # pulseseriesList=['OfflineIceTopHLCTankPulses']
							)


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

def eventFilter(frame):
	"""simple event filter"""
	if abs(frame["I3EventHeader"].event_id) % 100 > 10:
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
# tray.AddModule(eventFilter,"evtFilter",
# 				streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
#  				)
tray.AddModule(calcWeight,"calcWeight",
				streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
				)
tray.AddModule(icetop_Level3_scripts.modules.CheckContainment,
							 'MCPrimaryCheckContainment', Particles=["MCPrimary"], Detector="IC86.2019")


tray.AddModule("I3Writer","i3writer",
						# filename=str(outputDir)+"/hadronTimeTest/"+str(fileName)+"HadronNonTrigImProperTEvts.i3.gz",
						# filename=str(outputDir)+"/dataSetCleanTestSeedSame/"+str(fileName)+"CleanVEMEvts.i3.gz",
						filename=str(outputDir)+dataSetClean+str(fileName)+"CleanVEMEvts.i3.gz",
						# filename=str(outputDir)+"/dataSetCleanOfficial/"+str(fileName)+"CleanVEMEvts.i3.gz",
						streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						# streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
						# streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
						)

# tray.AddModule("I3NullSplitter","nullsplitter")

tray.Add(hdfwriter.I3HDFWriter, 'hdfNull',
		Output=str(outputDir)+dataSetClean+str(fileName)+"CleanVEMEvts.hdf5",
		# Output=str(outputDir)+"/dataSetCleanTestSeedSame/"+str(fileName)+"CleanVEMEvts.hdf5",
		# Output=str(outputDir)+"/dataSetCleanOfficial/"+str(fileName)+"CleanVEMEvts.hdf5",
		# Output=str(outputDir)+"/hadronTimeTest/"+str(fileName)+"NullHadronNonTrigImProperTEvts.hdf5",
		CompressionLevel=9,
		# SubEventStreams=['IceTopSplit'],
		SubEventStreams=['NullSplit'],
		# SubEventStreams=["nullsplitter",'IceTopSplit',"nullsplitter",'NullSplit',]
		# SubEventStreams=["ice_top"],
		# Streams=[icetray.I3Frame.DAQ],
		keys = [
		"MCPrimary","I3EventHeader","HLC6_5000","tank6_5000","tank6_4000","tank6_3000","tank6_2000",
		"tank7_5000","tank7_4000","tank7_3000","tank7_2000","tank8_5000","tank8_4000","tank8_3000","tank8_2000",
		"tank9_5000","tank9_4000","tank9_3000","tank9_2000","tank10_5000","tank10_4000","tank10_3000","tank10_2000",
		"OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses",
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
		"IceTopSTA5_13_filter","SDST_IceTopSTA3_13_filter","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank6",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank7","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank8",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_6SLCTankHit_t","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_7SLCTankHit_t",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_8SLCTankHit_t","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_6SLCTankHit_tSC",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_7SLCTankHit_tSC","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_8SLCTankHit_tSC",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank5","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank4"
		,"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_isSLCTank3","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_5SLCTankHit_t",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_4SLCTankHit_t","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_3SLCTankHit_t",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_2SLCTankHit_t","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_2SLCTankHit_tSC",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_5SLCTankHit_tSC","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_4SLCTankHit_tSC",
		"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_3SLCTankHit_tSC","HLCSLC_2TankHit_t","HLCSLC_2TankHit_tSC","HLCSLC_3TankHit_t",
		"HLCSLC_3TankHit_tSC","HLCSLC_4TankHit_t","HLCSLC_4TankHit_tSC","HLCSLC_5TankHit_t","HLCSLC_5TankHit_tSC","HLCSLC_6TankHit_t",
		"HLCSLC_7TankHit_tSC","HLCSLC_8TankHit_t","HLCSLC_6TankHit_tSC","HLCSLC_7TankHit_t","HLCSLC_8TankHit_tSC",
		"HLCSLC_hitStations","HLCSLC_hitTanks","HLCSLC_isTank1","HLCSLC_isTank2","HLCSLC_isTank3","HLCSLC_isTank4","HLCSLC_isTank5","HLCSLC_isTank6",
		"HLCSLC_isTank7","HLCSLC_isTank8","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_NoBackground","OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_NoBackground",
		"OfflineIceTopHLCVEMPulses_isSTA2","OfflineIceTopHLCVEMPulses_isSTA3","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_NoBackground_isSTA1",
		"OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_NoBackground_isSTA2","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_NoBackground_isSTA3",
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