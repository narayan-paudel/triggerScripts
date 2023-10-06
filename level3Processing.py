#!/usr/bin/env python3

'''
This example shows how to include IceTop DOMSets to the GCD file in addition to default InIce DOMSets
'''
import os
import sys

from icecube.icetray import I3Tray, I3Units
from icecube import dataclasses,icetray,dataio
from icecube.trigger_sim import GetDefaultDOMSets
from icecube.trigger_sim.InjectDefaultDOMSets import InjectDefaultDOMSets
# from icecube.filterscripts.offlineL2 import LaputopStandard
# from icecube.toprec import LaputopStandard
# from icecube.toprec import LaputopSmallShower
from icecube.icetop_Level3_scripts.segments.level3_IceTop import level3_IceTop
from icecube.recclasses import I3LaputopParams, LaputopParameter as Par

# from icecube.filterscripts.offlineL2 import Globals
# from icecube.filterscripts.offlineL2.level2_Reconstruction_Muon import OfflineMuonReco
# from icecube.filterscripts.offlineL2.level2_Reconstruction_IceTop import ReconstructIceTop
# from icecube.phys_services.which_split import which_split



import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec

from icecube import icetop_Level3_scripts
import icecube.icetray
from icecube import phys_services, sim_services
from icecube import tableio, hdfwriter
from icecube import topeventcleaning, tpx,toprec
# from icecube import top_background_simulator

import numpy as np
import math

from weighting.python.fluxes import GaisserH4a_IT
# from weighting.fluxes import GaisserH4a_IT
from weighting.python.weighting import icetop_mc_weights

weight_file = "/home/enpaudel/icecube/triggerStudy/simFiles/dataset_info200UnevenSnow.json" #my sim json file

ConfigIDMap = {102:(6,5000,"HLC6"),30043:(7,3000,"HG7")}
#config id mapped with threshold,timewindow and label
print("keys",ConfigIDMap.keys())

CORSIKA_ID = "DAT059871"
# outputDir = "/home/enpaudel/icecube/triggerStudy/simFiles/"
outputDir = "/data/sim/IceTop/2023/generated/untriggered/level3/"


# GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"
GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HGDomsets.i3.gz"
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
primary = args.input[0].split("/")[-1][0]

#example python level3Processing.py /home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique1_6/pDAT000027GenDetFiltProcUnique.i3.bz2

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

fileName = args.input[0].split("/")[-1]
fileName = fileName.split(".")[0]
fileName = fileName.split("Gen")[0]
print(fileName)
# print("inputs",args.input)

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
    # print("total no of events",self.totalEvents)
    for ikey in ConfigIDMap.keys():
      frame = self.checkTriggers(frame,triggerHierarchy, ikey)
    self.PushFrame(frame)

  def Finish(self):
    print("total number of events",self.totalEvents)

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
  # print("H4aWeight,type,zen energy flux,gen",p_weights,p_type,p_zenith,p_energy,flux(p_energy,p_type),generator(p_energy,p_type,np.cos(p_zenith)))

name = ""

tray = I3Tray()
tray.AddModule("I3Reader","reader",
             filenameList=[GCD]+args.input,
            )

# tray.AddSegment(ReconstructIceTop, 'ReconstructIceTop',
#     Pulses      = 'CleanedHLCTankPulses',
#     # Pulses      = 'IceTopTankPulses',
#     CoincPulses = 'CleanedCoincOfflinePulses',
#     If = which_split(split_name=Globals.filter_globals.IceTopSplitter)#'ice_top')
# )

tray.AddModule(TriggerCheck, "trigChk",
        # Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
        )
tray.AddModule('I3TankPulseMerger',
             # filenameList=[GCD]+args.input,
             InputVEMPulses = 'OfflineIceTopSLCVEMPulses',
             OutputTankPulses = 'OfflineIceTopSLCTankPulses',
             ExcludedTanks  = 'ExcludedSLCTanks'
             )
tray.AddModule(calcWeight,"calcWeight",
        streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics]
        )

tray.AddSegment(level3_IceTop,"level3",
  detector="IC86",
  do_select=False,
  isMC=True,
  simulate_background_rate=0,
  add_jitter=False,
  snowLambda=2.6,
  icetopStreamName="IceTopSplit",
  nameChanges=[["QTriggerHierarchy","I3TriggerHierarchy"]],
  ignore_slc_calib=True
  )

tray.AddModule("I3Writer","i3writer",
            filename=str(outputDir)+str(fileName)+"L3.i3.gz",
            streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            )

tray.Execute()
tray.Finish()
