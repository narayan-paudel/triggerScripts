#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import glob

from icecube import dataio,dataclasses, icetray, topeventcleaning, toprec
from I3Tray import I3Tray

inputPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/"

inputFiles = sorted(glob.glob(inputPath+"/*GenDetFiltProcUniqueCleanVEMEvts.i3.gz"))[:1]
print("input files",inputFiles)
outputFolder = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetInclFilt/"

GCD = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"

def inclinationFilter(frame,zenithThreshold):  
  if frame.Stop == icetray.I3Frame.Physics:
    if frame["I3EventHeader"].sub_event_stream in ['IceTopSplit']:
      zenith_true = frame["MCPrimary"].dir.zenith
      zenith_reco = frame["ShowerPlane"].dir.zenith
      if not np.isnan(zenith_reco):
        zenith_diff = zenith_reco - zenith_true
        if zenith_diff*180.0/np.pi >= zenithThreshold:
          frame["inclinedFilter"] = dataclasses.I3Double(1)
        elif zenith_diff*180.0/np.pi >= zenithThreshold:
          frame["inclinedFilter"] = dataclasses.I3Double(0)
      else:
        frame["inclinedFilter"] = dataclasses.I3Double(-1000)




tray = I3Tray()

tray.Add('I3Reader',
  FilenameList = [GCD] + inputFiles
  )

tray.Add('Delete',
  Keys = ["ShowerCOG","ShowerPlane","ShowerPlaneParams",'ClusterCleaningExcludedTanks',"ClusterCleaningExcludedTanksSLC"]
  )

def Unify(frame, Keys, Output):
  """
  Simple utility to merge RecoPulseSerieses into a single Union.
  """
  extants = [k for k in Keys if k in frame]
  union = dataclasses.I3RecoPulseSeriesMapUnion(frame, extants)
  frame[Output] = union

tray.Add(Unify,"UnionHLCSLC",
  Keys=["OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulses"],
  Output='IceTopTankPulses'
  )

tray.Add('I3TopRecoCore',
  DataReadout     = 'IceTopTankPulses',
  NTanks          = 7,
  ShowerCore      = 'ShowerCOG',
  Verbose         = False,
  Weighting_Power = 0.5,
  If              = lambda frame: 'IceTopTankPulses' in frame
  )

tray.Add('I3TopRecoPlane',
  DataReadout = 'IceTopTankPulses',
  ShowerPlane = 'ShowerPlane',
  Trigger     = 3,
  Verbose     = False,
  If          = lambda frame: 'IceTopTankPulses' in frame
)

# tray.Add(inclinationFilter,"inclFilt",
#   zenithThreshold = 20)

keep_all = ["H4aWeight","HLC6_5000","tank7_3000","I3EventHeader","OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulses",
'IceTopTankPulses',"ShowerCOG","ShowerPlane","ShowerPlaneParams","MCPrimary"]
# def keep(frame):

tray.AddModule( "Keep", "CleanUpKeys",
   Keys = keep_all)



tray.Add('I3Writer',
  Filename = outputFolder+"showerInclinedFilter.i3.gz",
  Streams  = [icetray.I3Frame.DAQ,
        icetray.I3Frame.Physics]
  )
tray.Add(hdfwriter.I3HDFWriter, 'hdfNull',
    Output=outputFolder+"showerInclinedFilter.hdf5",
    CompressionLevel=9,
    # SubEventStreams=['IceTopSplit'],
    SubEventStreams=['NullSplit'],
    # SubEventStreams=["nullsplitter",'IceTopSplit',"nullsplitter",'NullSplit',]
    # SubEventStreams=["ice_top"],
    # Streams=[icetray.I3Frame.DAQ],
    keys = ["H4aWeight","HLC6_5000","tank7_3000","I3EventHeader","OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulses",
'IceTopTankPulses',"ShowerCOG","ShowerPlane","ShowerPlaneParams","MCPrimary"
    ]
    )

tray.Execute()
tray.Finish()
