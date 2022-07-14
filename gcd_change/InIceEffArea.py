#!/usr/bin/env python

#This calculates the effective area of the in-ice volume as a function of zenith angle


from icecube import radcube, dataclasses, dataio, icetray
from icecube.icetray import I3Units

import numpy as np
from scipy.spatial import ConvexHull
import os
ABS_PATH_HERE=str(os.path.dirname(os.path.realpath(__file__)))

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--gcd', type=str, default=ABS_PATH_HERE+"/data/GCDMC-ITDec30Gen2ScintAntennasLifted_OutriggerFullCirle.i3.gz", help='GCD file for the event')
parser.add_argument('--gen1', action="store_true", help='If set, will not consider scints')
args = parser.parse_args()

def ParseGCDFile(gcdFilename):
  #Reads in the Gframe, parses the in-ice locations
  #Returns the location of the 3d verticies that define the hull

  print("Parsing GCD file", gcdFilename)
  file = dataio.I3File(gcdFilename, 'r')
  frame = file.pop_frame()
  while frame.Stop != icetray.I3Frame.Geometry and file.more():
    frame = file.pop_frame()

  assert(frame.Stop == icetray.I3Frame.Geometry) #GCD file must have GFrame
  assert(frame.Has("I3Geometry")) #Needs to have in-ice info
  geometry = frame["I3Geometry"]
  omgeomap = geometry.omgeo

  omList = []
  for key in omgeomap.keys():
    if key.om in [61, 62, 63, 64]: #Skip IT tanks
      continue
    if args.gen1 and key.string > 86:
      continue

    om = omgeomap[key]
    omList.append([om.position.x, om.position.y, om.position.z])


  print("Read in {} in-ice OMs".format(len(omList)))

  assert(len(omList) >= 4) #Must have enough in-ice doms to do a 3D hull

  omList = np.array(omList)

  hull3d = ConvexHull(omList)
  mesh3d = omList[hull3d.vertices] #Get the locations defining the hull

  return mesh3d


def GetProjectionArea(mesh3d, direction):
  #Given the location of the 3d verticies that define the hull and the direction of the shower
  #Finds the projection of the volume for the specificed arrival direction
  #Returns the projected area

  #Convert to I3Position and get the hull in shower plane coords
  i3Mesh = [dataclasses.I3Position(pos[0], pos[1], pos[2]) for pos in mesh3d]
  i3MeshShower = [radcube.GetShowerFromIC(pos, direction) for pos in i3Mesh]

  meshShower = np.array([[pos.x, pos.y] for pos in i3MeshShower]) #Compress to 2D
  hull2d = ConvexHull(meshShower)
  mesh2d = meshShower[hull2d.vertices]


  lines = np.hstack([mesh2d,np.roll(mesh2d,-1,axis=0)])
  area = 0.5*abs(sum(x1*y2-x2*y1 for x1,y1,x2,y2 in lines))
  return area


mesh3d = ParseGCDFile(args.gcd)


height = max(mesh3d[:,2]) - min(mesh3d[:,2])
areaVertical = GetProjectionArea(mesh3d, dataclasses.I3Direction(180, 0))

print("Height is", height / I3Units.km, " km")
print("Area is", areaVertical / I3Units.km2, "km^2")
print("Volume is", height * areaVertical / (I3Units.km)**3, " km^3")


cosZen = [0.05, 0.10, 0.15, 0.20, 0.25, 0.3]

print("Cos(theta)\tAeff [km^2]")

for cos in cosZen:
  area = 0
  nSteps = 100

  #Average over phi (even though it's pretty symmetric)
  for phi in np.linspace(0, 2*np.pi, nSteps):
    direction = dataclasses.I3Direction(np.arccos(cos), phi)
    area += GetProjectionArea(mesh3d, direction)

  area /= nSteps

  print(cos, "\t\t{0:0.4f}".format(area / icetray.I3Units.km2))