#!/bin/env python3

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio

from icecube import tableio, hdfwriter
from icecube import topeventcleaning, tpx,toprec

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default=
  "/data/sim/IceCubeUpgrade/CosmicRay/Narayan/dataSetCleanTanks/pDAT005988GenDetFiltProcUniqueCleanVEMEvts.i3.gz",
   help='Input files after running detector.py.')
args = parser.parse_args()

GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"

print(args.input)

# outFile = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetFullEff/"+args.input[0].split("/")[-1]
outFile = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetFullEffBoth/"+args.input[0].split("/")[-1]
print(outFile)

relevant_keeps = ["tank7_3000","SMT102","SMT273","QTriggerHierarchy","QFilterMask",
"OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge",
"OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopHLCTankPulsesCleanTimeCleanCharge",
"MCPrimary","IceTopRawData","IceTopPulses","I3Triggers","I3EventHeader","H4aWeight","HLC6_5000"
]

#example python selectFullEfficiency.py 
#/data/sim/IceCubeUpgrade/CosmicRay/Narayan/dataSetCleanTanks/pDAT005988GenDetFiltProcUniqueCleanVEMEvts.i3.gz

def energySelect(frame):
  p_energy = frame["MCPrimary"].energy*I3Units.GeV/I3Units.eV
  p_type = frame["MCPrimary"].type
  p_zenith = frame["MCPrimary"].dir.zenith
  # print("p_energy",p_energy,I3Units.GeV,I3Units.eV)
  if 10**16 <= p_energy <= 10**16.6:
    # print("p_energy",p_energy,I3Units.GeV,I3Units.eV)
    return True
  else:
    return False

def triggerSelect(frame):
  if frame.Has("SMT273"):
    # if frame["SMT273"] > icetray.I3Int(0) and frame["SMT102"] == icetray.I3Int(0):
    if frame["SMT273"] > icetray.I3Int(0) or frame["SMT102"] > icetray.I3Int(0):
      # print("SMT",frame["SMT273"],frame["SMT102"])
      return True
    else:
      return False
  else:
    return False

def subEventSelect(frame):
  if frame.Has("I3EventHeader"):
    if frame["I3EventHeader"].sub_event_stream == "convertedP":
      print(frame["I3EventHeader"].sub_event_stream)
      return True
    else:
      return False

tray = I3Tray()
tray.AddModule("I3Reader","reader",
  filenameList=[GCD]+args.input)

# tray.Add(energySelect,"selectE",
#   streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])
tray.Add(triggerSelect,"selectTrig",
  streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])
tray.Add(subEventSelect,"selectSubEvt",
  streams=[icetray.I3Frame.Physics])

tray.AddModule("Keep", "keep_relevant",
               keys = relevant_keeps,
               If=lambda x:True
               )
def Unify(frame, Keys, Output):
  """
  Simple utility to merge RecoPulseSerieses into a single Union.
  """
  extants = [k for k in Keys if k in frame]
  union = dataclasses.I3RecoPulseSeriesMapUnion(frame, extants)
  frame[Output] = union

tray.Add(Unify,"UnionHLCSLC",
  Keys=["OfflineIceTopHLCTankPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge"],
  Output='IceTopTankPulses',
  streams=[icetray.I3Frame.DAQ]
  )
tray.Add("I3NullSplitter",
  SubEventStreamName="convertedP")

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

tray.AddModule("I3Writer","i3writer",
            filename=outFile,
            streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            )
tray.Execute()
tray.Finish()