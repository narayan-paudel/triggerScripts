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





tray = I3Tray()
tray.AddModule("I3Reader","reader",
	           filenameList=[GCD]+args.input,
	           # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
	           # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
	          )

# tray.AddModule("I3Writer","i3writer",
# 	          filename=str(outputDir)+"trigChkTestSeparate/"+str(fileName)+"TrigCount.i3.gz",
# 	          # streams=[icetray.I3Frame.DAQ],
# 	          streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
# 	          # streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus]
# 	          )

# tray.AddModule("I3NullSplitter","nullsplitter")

tray.Add(hdfwriter.I3HDFWriter, 'hdfNull',
    Output=str(outputDir)+"/trigChkTestSeparate/"+str(fileName)+"TrigCountNullSplitSingle.hdf5",
    CompressionLevel=9,
    # SubEventStreams=['IceTopSplit'],
    SubEventStreams=['NullSplit'],
    # SubEventStreams=["nullsplitter",'IceTopSplit',"nullsplitter",'NullSplit',]
    # SubEventStreams=["ice_top"],
    # Streams=[icetray.I3Frame.DAQ],
    keys = [
    "MCPrimary","ITSMTTriggered","ITGlobalTriggered","OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulsesTotalCharge",
    "OfflineIceTopHLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalHit","OfflineIceTopSLCTankPulsesTotalHit","OfflineIceTopHLCVEMPulses",
    "OfflineIceTopSLCVEMPulses","OfflineIceTopSLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalHit",
    "OfflineIceTopSLCVEMPulsesTotalHit","OfflineIceTopSLCTankPulsesHitTimeDuration","OfflineIceTopHLCTankPulsesHitTimeDuration",'ITCleanSLCTankPulsesHitTimeDuration','ITCleanSLCTankPulsesTotalHit'
    ]
    )

tray.Add(hdfwriter.I3HDFWriter, 'hdfIce',
    Output=str(outputDir)+"/trigChkTestSeparate/"+str(fileName)+"TrigCountIceTopSplitSingle.hdf5",
    CompressionLevel=9,
    SubEventStreams=['IceTopSplit'],
    # SubEventStreams=['NullSplit'],
    # SubEventStreams=["nullsplitter",'IceTopSplit',"nullsplitter",'NullSplit',]
    # SubEventStreams=["ice_top"],
    # Streams=[icetray.I3Frame.DAQ],
    keys = [
    "MCPrimary","ITSMTTriggered","ITGlobalTriggered","OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses","OfflineIceTopSLCTankPulsesTotalCharge",
    "OfflineIceTopHLCTankPulsesTotalCharge","OfflineIceTopHLCTankPulsesTotalHit","OfflineIceTopSLCTankPulsesTotalHit","OfflineIceTopHLCVEMPulses",
    "OfflineIceTopSLCVEMPulses","OfflineIceTopSLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalCharge","OfflineIceTopHLCVEMPulsesTotalHit",
    "OfflineIceTopSLCVEMPulsesTotalHit","OfflineIceTopSLCTankPulsesHitTimeDuration","OfflineIceTopHLCTankPulsesHitTimeDuration",'ITCleanSLCTankPulsesHitTimeDuration','ITCleanSLCTankPulsesTotalHit'
    ]
    )

tray.Execute()
tray.Finish()
