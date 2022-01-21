#!/usr/bin/python

import sys
from I3Tray import *
from icecube import icetray, dataclasses,recclasses, phys_services,dataio, rootwriter, hdfwriter, tableio, tpx

from icecube.weighting.fluxes import GaisserH4a_IT
from icecube.weighting.weighting import icetop_mc_weights

from os.path import expandvars
import numpy as np

rootout = sys.argv[1]
infiles = sys.argv[2:]

tray = I3Tray()

tray.AddModule('I3Reader', 'reader',FileNameList = infiles)

flux = GaisserH4a_IT()
'''
#LE = icetop_mc_weights(12360, dataset_file="/data/ana/CosmicRay/MuonDensity/level3_mc/output/dataset_info_new.json") # H
#HE = icetop_mc_weights(20143, dataset_file="/data/ana/CosmicRay/MuonDensity/level3_mc/output/dataset_info_new.json") # H
#LE = icetop_mc_weights(12630, dataset_file="/data/ana/CosmicRay/MuonDensity/level3_mc/output/dataset_info_new.json") # He
#HE = icetop_mc_weights(20145, dataset_file="/data/ana/CosmicRay/MuonDensity/level3_mc/output/dataset_info_new.json") # He
#LE = icetop_mc_weights(12631, dataset_file="/data/ana/CosmicRay/MuonDensity/level3_mc/output/dataset_info_new.json") # O
#HE = icetop_mc_weights(20146, dataset_file="/data/ana/CosmicRay/MuonDensity/level3_mc/output/dataset_info_new.json") # O
#LE = icetop_mc_weights(12362, dataset_file="/data/ana/CosmicRay/MuonDensity/level3_mc/output/dataset_info_new.json") # Fe
#HE = icetop_mc_weights(20144, dataset_file="/data/ana/CosmicRay/MuonDensity/level3_mc/output/dataset_info_new.json") # Fe
'''

#LE = icetop_mc_weights(12360, dataset_file="/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json") #H
#HE = icetop_mc_weights(20143, dataset_file="/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json") #H
#LE = icetop_mc_weights(12630, dataset_file="/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json") #He
#HE = icetop_mc_weights(20145, dataset_file="/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json") #He
#LE = icetop_mc_weights(12631, dataset_file="/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json") #O
#HE = icetop_mc_weights(20146, dataset_file="/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json") #O
LE = icetop_mc_weights(12362, dataset_file="/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json") #Fe
HE = icetop_mc_weights(20144, dataset_file="/data/ana/CosmicRay/IceTop_level3/sim/L3_MC_incomplete_test/output/dataset_info.json") #Fe

generator = LE + HE

def mcprim(frame):
    p_energy = frame["MCPrimary"].energy 
    p_type = frame["MCPrimary"].type
    p_zenith = frame["MCPrimary"].dir.zenith
    p_weights = flux(p_energy,p_type)/generator(p_energy,p_type,np.cos(p_zenith))
    frame["H4awei"]=dataclasses.I3Double(p_weights)
   # print p_weights
tray.AddModule(mcprim,'mcprim')

def itstdfilt(frame):
    return ('IT73AnalysisIceTopQualityCuts' in frame and all(frame["IT73AnalysisIceTopQualityCuts"].values()))
tray.AddModule(itstdfilt,'itstdfilt')


#hdf = hdfwriter.I3HDFTableService(hdfout)
root = rootwriter.I3ROOTTableService(rootout, 'MasterTree')
tray.AddModule(tableio.I3TableWriter, 'TableWriter',
#    TableService    =[hdf],
    TableService    =[root],
    SubEventStreams = ['ice_top'],
    Keys            = [
                       #'IceTopHLCPulseInfo','CleanedHLCTankPulses',
        'MCPrimary','H4awei','I3EventHeader','IT73AnalysisIceTopQualityCuts',
        'QFilterMask','Laputop','LaputopParams','LaputopSmall','LaputopSmallParams',
        'OfflineIceTopHLCTankPulses','OfflineIceTopSLCTankPulses',
        'IceTopHLCSeedRTPulses',
        'IceTopLaputopSeededSelectedHLC','IceTopLaputopSeededSelectedSLC',
        dict(key = 'OfflineIceTopHLCTankPulses', converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(), name = 'HLCTankPulses_EventInfo'),
            dict(key = 'OfflineIceTopSLCTankPulses', converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(), name = 'SLCTankPulses_EventInfo'),

             dict(key = 'IceTopHLCSeedRTPulses', converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(), name = 'SeedRTHLC_EventInfo'),

             dict(key = 'IceTopLaputopSeededSelectedHLC', converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(), name = 'LapuSeedHLC_EventInfo'),
            dict(key = 'IceTopLaputopSeededSelectedSLC', converter = phys_services.converters.I3EventInfoConverterFromRecoPulses(), name = 'LapuSeedSLC_EventInfo')
                       ] )
tray.Execute()
tray.Finish()

