#!/usr/bin/env python3

'''
This example shows how to include IceTop DOMSets to the GCD file in addition to default InIce DOMSets
'''
import os
import sys

from icecube.icetray import I3Tray, I3Units
from icecube import dataclasses,icetray,dataio,simclasses

from glob import glob

import simweights

import numpy as np

from weighting.python.fluxes import GaisserH4a_IT
# from weighting.fluxes import GaisserH4a_IT
from weighting.python.weighting import icetop_mc_weights

# import matplotlib.pyplot as plt
import pylab as plt

weight_file = "/home/enpaudel/icecube/triggerStudy/simFiles/dataset_info200UnevenSnow.json" #my sim json file


# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/ITGenTest/FeDAT000001GenDet.i3.bz2", help='Input files after running detector.py.')
# args = parser.parse_args()

FILE_DIR = "/home/enpaudel/icecube/triggerStudy/simFiles/level3/level3"
filelist = sorted(glob(FILE_DIR + "/pDAT00000*.i3.gz"))
assert filelist

# create a dictionary that mimics the structure of a pandas/h5py table
I3TopInjectorInfo = {k: [] for k in dir(simclasses.I3TopInjectorInfo) if k[0] != "_"}

# Same for I3CorsikaWeight but we only need a few columns
I3CorsikaWeight: dict = {"energy": [], "type": [], "zenith": [], "weight": []}
MCPrimary: dict = {"energy": [], "type": [], "zenith": []}

# MCPrimary = {k: [] for k in dir(dataclasses.I3Particle) if k[0] != "_"}
print(MCPrimary)

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

conv_weights = []


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
  # frame["H4aWeight"]=dataclasses.I3Double(p_weights)
  # print("H4aWeight,type,zen energy flux,gen",p_weigh
  return p_weights



# loop over all the files we want to read
for filename in filelist:
    # open the i3 files with the dataio interface
    infile = dataio.I3File(filename)
    print("Reading " + filename)

    # loop over the frames in the file
    while infile.more():
        # get the frame
        frame = infile.pop_frame()


        # if this is an S-Frame
        if frame.Stop == frame.Simulation:
            # get the info from the frame
            info = frame["I3TopInjectorInfo"]

            for k in I3TopInjectorInfo:
                I3TopInjectorInfo[k].append(getattr(info, k))

        # if this is a physics event in the right sub-event stream
        elif frame.Stop == frame.Physics and frame["I3EventHeader"].sub_event_stream == "IceTopSplit":
          conv_weight = calcWeight(frame)
          conv_weights.append(conv_weight)
          # get the weighting object
          w = frame["MCPrimary"]
          # for each of the columns we need get it from the frame object
          # and put it in the correct column
          MCPrimary["energy"].append(w.energy)
          MCPrimary["type"].append(w.type)
          MCPrimary["zenith"].append(w.dir.zenith)
          # MCPrimary["weight"].append(w.weight)

# make a dictionary object to mimic the file structure of a pandas file
fileobj = {"I3TopInjectorInfo": I3TopInjectorInfo,"MCPrimary":MCPrimary}

# create the weighter object
weighter = simweights.IceTopWeighter(fileobj)

# create an object to represent our cosmic-ray primary flux model
flux = simweights.GaisserH4a_IT()

# get the weights by passing the flux to the weighter
simWeightList = weighter.get_weights(flux)



# print some info about the weighting object
# print(weighter.tostring(flux))

# create equal spaced bins in log space
# bins = plt.geomspace(3e4, 1e6, 50)

# get energy of the primary cosmic-ray from `PolyplopiaPrimary`
primary_energy = weighter.get_weight_column("energy")
# print("primary_energy",primary_energy)


print("conv_weight",len(conv_weights),conv_weights[:2])
print("sim weights",len(simWeightList),simWeightList[:2])
# histogram the primary energy with the weights
# plt.hist(simWeightList)
bins = plt.geomspace(3e4, 1e6, 50)
plt.hist(primary_energy, weights=simWeightList, bins=bins,histtype="step",label="simweights")
plt.hist(primary_energy, weights=conv_weights, bins=bins,histtype="step",label="weighting")
# make the plot look good
plt.xlabel("Primary Energy [GeV]")
plt.ylabel("Event Rate [Hz]")
plt.yscale("log")
# plt.xlim(bins[0], bins[-1])
# plt.ylim(0.1, 10)
# plt.loglog()
plt.legend()
plt.savefig(plotFolder+"/weighting_compare.pdf")
plt.show()

