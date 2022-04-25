#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

from scipy.spatial import ConvexHull

from icecube import dataclasses, icetray, dataio, radcube, rock_bottom
from I3Tray import *

import numpy as np

import os
ABS_PATH_HERE=str(os.path.dirname(os.path.realpath(__file__)))

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Input data files.')
parser.add_argument('--output', type=str, default=ABS_PATH_HERE+"/../plots/GCDPlot.pdf", help='Output data files.')
parser.add_argument('--phase1', action="store_true", help='Only use the Gen1 scints/DOMs.')
parser.add_argument('--phase2', action="store_true", help='Only use the Gen1 scints/DOMs.')
parser.add_argument('--showids', action="store_true", help='Show IDs')

args = parser.parse_args()

framenames=["I3AntennaGeometry", "I3Geometry", "I3ScintGeometry", "I3IceActGeometry"]

def AddAntennas(array, i3geometry, numLabel=None):
  antennageomap = i3geometry.antennageo
  for key in antennageomap.keys():

    if args.phase1 and key.station > 10032:
      continue

    if args.phase2 and not (key.station <= 10032 or key.station in phase2Stations):
      continue

    pos = antennageomap[key].position
    array.append(pos)

def AddIceAct(array, i3geometry, numLabel=None):
  iceactgeomap = i3geometry.iceactgeo
  for key in iceactgeomap.keys():

    pos = iceactgeomap[key].position
    array.append(pos)

phase2Stations = [10033, 10034, 10035, 10036, 10037, 10038, 10039, 10041, 10042, 
                  10043, 10044, 10048, 10049, 10050, 10055, 10056, 10057, 10058, 
                  10065, 10066, 10067, 10068, 10076, 10077, 10078,
                  10155, 10156, 10157, 10158, 10159, 10160] #This line is the infill stations

def AddScints(array, i3geometry, numLabel=None):
  scintgeomap = i3geometry.scintgeo
  for key in scintgeomap.keys():

    if args.phase1 and key.station > 10032:
      continue

    if args.phase2 and not (key.station <= 10032 or key.station in phase2Stations):
      continue

    pos = scintgeomap[key].position
    array.append(pos)
    if numLabel != None and not key.station in numLabel.keys():
      numLabel[key.station] = [pos.x, pos.y]

def AddStns(array, snow, i3geometry, numLabel=None):
  stationgeomap = i3geometry.stationgeo
  for key in stationgeomap.keys():
    stn = stationgeomap[key]
    for tank in stn:
      pos = tank.position
      array.append(pos)
      snow.append(dataclasses.I3Position(pos.x, pos.y, tank.snowheight+pos.z))

phase2Strings = [1001, 1003, 1005, 1007, 1009, 1010, 1012, 1015, 1017, 1020, 
                 1022, 1023, 1025, 1028, 1030, 1033, 1035, 1036, 1038, 1041, 
                 1043, 1046, 1118, 1119, 1120]

def AddInIce(array, i3geometry, numLabel=None):
  omgeomap = i3geometry.omgeo
  for key in omgeomap.keys():
    if key.om in [61, 62, 63, 64] and key.string <= 83: #Skip Icetop doms
      continue

    if args.phase1 and key.string > 83:
      continue

    if args.phase2 and not (key.string <= 83 or key.string in phase2Strings):
      continue

    om = omgeomap[key]
    
    # Ensure nothing on the surface
    if om.position.z > 1900:
      continue

    array.append([om.position.x, om.position.y, om.position.z])

    if numLabel != None and not key.string in numLabel.keys():
      numLabel[key.string] = [om.position.x, om.position.y]

def PlotList(ax, array, label, marker, color, alpha=1, hull=False, d3 = False, surf=False, numLabel=[], s= 0.7, **kargs):
  if not len(array):
    return

  x = []
  y = []
  z = []
  for pos in array:
    x.append(pos[0])
    y.append(pos[1])
    z.append(pos[2])

  if d3:
    print("{0} \tCenter({1:0.1f},{2:0.1f}) Max:{3:0.1f}, Min:{4:0.1f}".format(label, np.average(x), np.average(y), max(z), min(z)))


  if surf:
    ax.plot_trisurf(x, y, z, label=label, color=color, alpha = 0.3, shade=False)
  elif d3:
    ax.scatter(x, y, z, label=label, marker=marker, alpha=alpha, color=color, s=s, **kargs)
  else:
    ax.scatter(x, y, label=label, marker=marker, alpha=alpha, color=color, s=s, **kargs)

  if len(numLabel) and args.showids:
    print("Plotting labels")
    for key in numLabel.keys():
      pos = numLabel[key]
      ax.text(pos[0], pos[1], str(key))


file = dataio.I3File(args.input, "r")
gFrame = file.pop_frame()
gFrame = file.pop_frame()

antList = []
iceactList = []
scintList = []
scintNums = {}
stnList = []
stnNums = {}
snow = []
iniceList = []
iniceNums = {}

print("Reading in data from", args.input)
for name in framenames:
  if not name in gFrame:
    print("There is no", name, "in the GFrame")
    print(gFrame)
    continue

  AddAntennas(antList, gFrame[name])
  AddScints(scintList, gFrame[name], scintNums)
  AddStns(stnList, snow, gFrame[name])
  AddInIce(iniceList, gFrame[name], iniceNums)
  AddIceAct(iceactList, gFrame[name])

print("N ants:",len(antList), "\nN scints:", len(scintList), "\nN tanks:", len(stnList), "\nN IceAct:", len(iceactList))

fig =plt.figure(figsize=(18, 8))
gs = gridspec.GridSpec(1, 2, wspace=0.0, hspace=0.0, width_ratios=[1,2])
ax = fig.add_subplot(gs[0], projection='3d')

PlotList(ax, stnList, "IT", 'x', 'k', d3=True)
PlotList(ax, snow, "Snow", 's', 'k', d3=True, surf=True)
PlotList(ax, scintList, "Scint.", 'o', 'b', d3=True)
PlotList(ax, antList, "Ant.", 's', 'r', d3=True)
PlotList(ax, iceactList, "IceAct", '^', 'g', d3=True, surf=True)
# PlotList(ax, iniceList, "Inice", 's', 'k', d3=True, alpha = 0.2)

############################

ax.set_aspect('equal', 'box')
ax.view_init(25,135-180)
ax.set_xlabel("East / m")
ax.set_ylabel("North / m")

ax = fig.add_subplot(gs[1])

PlotList(ax, stnList, "IT", 'x', 'k')
PlotList(ax, scintList, "Scint.", 'o', 'b', numLabel=scintNums)
PlotList(ax, antList, "Ant.", 's', 'r')
PlotList(ax, iceactList, "IceAct", '^', 'g', s=6)

ax.set_aspect('equal')
ax.set_xlabel("East / m")
ax.set_ylabel("North / m")

print("Making file", args.output)
plt.savefig(args.output, bbox_inches='tight')
plt.close()