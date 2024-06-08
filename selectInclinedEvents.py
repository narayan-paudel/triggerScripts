#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio

import numpy as np
import glob

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--vertEvents", dest="vertEvts",
                    default=False, action="store_true", required=False,
                    help="plot vertical events")
args = parser.parse_args()

#usage python selectInclinedEvents.py --vertEvents

inclinationCut = 60 #degree
energyCut = 10**16 #eV
fileDir = "/data/user/enpaudel/triggerStudy/simFiles/dataSetClean1_6/"

outputDirIncl = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/"
outputDirVert = "/home/enpaudel/dataExp/dataSetClean_VerticalLE/"

vertEvts = args.vertEvts

def testZenith(frame,zcut):
  return np.rad2deg(frame["MCPrimary"].dir.zenith) >= zcut

def testEnergy(frame,Ecut):
  return frame["MCPrimary"].energy * 10**9 >= Ecut

def selectEvents(frame,zcut,Ecut):
  return testZenith(frame,zcut) and testEnergy(frame,Ecut)

def selectVertEvents(frame,zcut,Emin,Emax):
  return np.rad2deg(frame["MCPrimary"].dir.zenith) <= zcut and frame["MCPrimary"].energy * 10**9 >= Emin and frame["MCPrimary"].energy * 10**9 < Emax

fileList = sorted(glob.glob(fileDir+"*.i3.*"))
for ifile in fileList:
  fileName = ifile.split("/")[-1]
  tray = I3Tray()
  tray.AddModule("I3Reader","reader",
               # filenameList = args.input,
               filename = ifile,
              )
  if not vertEvts:
    tray.AddModule(selectEvents,"selEvts",
                Ecut = 10**16,
                zcut = 60,
                streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

    tray.AddModule("I3Writer","i3writer",
                filename=str(outputDirIncl)+fileName,
                streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
                )
  else:
    tray.AddModule(selectVertEvents,"selVertEvts",
                Emax = 10**16,
                Emin = 10**15,
                zcut = 10,
                streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

    tray.AddModule("I3Writer","i3writer",
                filename=str(outputDirVert)+fileName,
                streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
                )


  tray.Execute()
  tray.Finish()