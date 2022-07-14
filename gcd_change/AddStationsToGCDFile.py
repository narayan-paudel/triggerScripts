#!/usr/bin/env python

#Reads in an existing GCD file with the optical strings and a second GCD file with the locations of the Gen1 stations
#Makes a new GCD file will the Gen1 stations, all the DOMs, and adds a station on top of each Gen2 string

#./AddStationsToGCDFile.py ../data/GCDMC-DOMOnly-Gen2_Sunflower_240m_v3.i3 --keepit

from icecube import icetray, dataclasses, dataio, radcube
from icecube.icetray import I3Units
import os
ABS_PATH_HERE=str(os.path.dirname(os.path.realpath(__file__)))
from os.path import expandvars
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Location of input GCD file to add stations to')
parser.add_argument('--output', type=str, default=ABS_PATH_HERE+"/../data/GCD-StationsAdded.i3.gz", help='Location of output GCD file after adding stations')
parser.add_argument('--gen1gcd', type=str, default=ABS_PATH_HERE+"/../data/GCD-Simulation-AntennaScint.i3.gz", help='GCD file to use for the Gen1 footprint')
parser.add_argument('--scintarmlength', type=float, default=72, help="The length of the arms from the center of the stations")
parser.add_argument('--scintpadding', type=float, default=5., help="Separation between pairs of panels")
parser.add_argument('--zval', type=float, default=2838.0 - dataclasses.I3Constants.OriginElev, help="The length of the arms from the center of the stations")
parser.add_argument('--spokeangles', nargs="+", default=[np.pi/2, 7.*np.pi/6., 11.*np.pi/6.], help="The angles of the spokes")
parser.add_argument('--keepit', action="store_true", help='If selected, IceTop will be kept in the file')
args = parser.parse_args()


def MakeScintGeo(name, position, orientation=0.0, heightAboveSnow=0.5):
  '''
  Makes a single instance of a scintillator geometry item
  Arguments:
  name - the name (str) of the tank (don't really matter)
  position - I3Position of the panel location
  orientation - angle in I3Units of the rotation angle of the panel
  heightAboveSnow - how high above the snow the panel is

  Returns:
  I3ScintGeo instance corresponding to the inputs
  '''
  scintGeo = dataclasses.I3ScintGeo()
  scintGeo.position = position
  up    = dataclasses.I3Direction(0,0,1)
  d_dir = dataclasses.I3Direction(1,0,0)
  d_dir.rotate_z(orientation*icetray.I3Units.rad)
  scintGeo.orientation = dataclasses.I3Orientation(d_dir, up)
  scintGeo.heightAboveSnow = heightAboveSnow
  scintGeo.scintName = name
  scintGeo.scintType = dataclasses.I3ScintGeo.ScintType.MADKIT

  return scintGeo

def MakeAntennaGeo(name, position, heightAboveSnow=0.5):
  '''
  Makes a single instance of a antenna geometry item
  Arguments:
  name - the name (str) of the tank (don't really matter)
  position - I3Position of the antenna location
  heightAboveSnow - how high above the snow the antennas is

  Returns:
  I3AntennaGeo instance corresponding to the inputs
  '''
  antennaGeo = dataclasses.I3AntennaGeo()
  antennaGeo.position = position
  antennaGeo.orientation = dataclasses.I3Orientation(1, 0, 0, 0, 0, 1)
  antennaGeo.heightAboveSnow = heightAboveSnow
  antennaGeo.antennaName = name
  antennaGeo.antennatype = dataclasses.I3AntennaGeo.AntennaType.SKALA2
  antennaGeo.cableLength = 50*I3Units.m

  return antennaGeo

def FindStringLocations(omgeomap):
  '''
  Finds the location of the non-Gen1 strings
  Arguments:
  omgeomap - An I3OMGeoMap containting all in-ice doms

  Returns:
  list of I3Positions corresponding to the non-Gen1 strings
  '''
  seenStrings = []
  stringLocations = []

  for key in omgeomap.keys():
    if key.string > 86 and not key.string in seenStrings:
      seenStrings.append(key.string)
      stringLocations.append(omgeomap[key].position)

  return stringLocations

def AddScintPositions(stationId, pos, scintgeomap):
  '''
  Given the center of a station, makes a N-spoked staition with 2 panels
  at the center and 2 at the end of each spoke. Spoke lengths are given by
  the arguments to this script as well as the spacing between neighboring panels

  Arguments:
  stationId - int corresponding to the station number (used for keys)
  pos - I3Position of the station center
  scintgeomap - I3ScintGeoMap that will be filled with I3ScintGeo instances
  '''

  #Add the stations in the center
  pos1 = dataclasses.I3Position(args.scintpadding/2.+pos.x, 0+pos.y, 0)
  pos1.z = args.zval
  pos2 = dataclasses.I3Position(-args.scintpadding/2.+pos.x, 0+pos.y, 0)
  pos2.z = args.zval
  scintgeomap[dataclasses.ScintKey(stationId, 1)] = MakeScintGeo("stn{}_panel1".format(stationId), pos1, 0.)
  scintgeomap[dataclasses.ScintKey(stationId, 2)] = MakeScintGeo("stn{}_panel2".format(stationId), pos2, 0.)

  #For each of the arm orientations, add one station
  for irot, rot in enumerate(args.spokeangles):
    pos1 = dataclasses.I3Position(args.scintarmlength, args.scintpadding/2., 0.)
    pos1.rotate_z(rot)
    pos1 += pos
    pos1.z = args.zval
    scintgeomap[dataclasses.ScintKey(stationId, 3+irot*2)] = MakeScintGeo("stn{}_panel{}".format(stationId, 3+irot*2), pos1, rot)

    pos2 = dataclasses.I3Position(args.scintarmlength, -args.scintpadding/2., 0.)
    pos2.rotate_z(rot)
    pos2 += pos
    pos2.z = args.zval
    scintgeomap[dataclasses.ScintKey(stationId, 4+irot*2)] = MakeScintGeo("stn{}_panel{}".format(stationId, 4+irot*2), pos2, rot)


def AddAntennaPositions(stationId, pos, antgeomap):
  '''
  Given the center of a station, makes a N-spoked staition with
  antennas half way down each spoke. Spoke lengths are given by
  the arguments to this script as well as the spacing between neighboring panels

  Arguments:
  stationId - int corresponding to the station number (used for keys)
  pos - I3Position of the station center
  antgeomap - I3AntennaGeoMap that will be filled with I3AntennaGeo instances
  '''
  for irot, rot in enumerate(args.spokeangles):
    posAnt = dataclasses.I3Position(args.scintarmlength/2., 0., 0.)
    posAnt.rotate_z(rot)
    posAnt += pos
    posAnt.z = args.zval
    antgeomap[dataclasses.AntennaKey(stationId, 3+irot*2)] = MakeAntennaGeo("stn{}_ant{}".format(stationId, 3+irot), posAnt)


def AddStation(stationId, pos, scintgeomap, antgeomap):
  AddScintPositions(stationId, pos, scintgeomap)
  AddAntennaPositions(stationId, pos, antgeomap)

def GetGen1Antennas():
  for frame in dataio.I3File(args.gen1gcd, "r"):
    if frame.Stop == icetray.I3Frame.Geometry:
      return frame[radcube.GetDefaultGeometryName()].antennageo

def GetGen1Scints():
  for frame in dataio.I3File(args.gen1gcd, "r"):
    if frame.Stop == icetray.I3Frame.Geometry:
      return frame["I3ScintGeometry"].scintgeo


def ParseGeometryFrame(gFrame):
  '''
  Reads in the info of the base GCD file, adds stations on top of each of the
  non-Gen1 strings. Adds all the stations from the Gen1 array. Adds all the
  DOMS back into the GFrame. Returns a new copy of the frame with updated info.

  Arguments:
  gFrame - Geometry frame from the baseline DOM layout

  Returns:
  Geometry frame with all the additional fields
  '''

  assert(frame.Stop == icetray.I3Frame.Geometry) #Make sure it is a G-Frame

  print("Finding strings...")
  stringLocations = FindStringLocations(frame["I3Geometry"].omgeo)
  print("There are Gen2 {} strings".format(len(stringLocations)))

  scintGeometry = dataclasses.I3Geometry()
  scintGeometry.scintgeo = dataclasses.I3ScintGeoMap(GetGen1Scints())
  antGeometry = dataclasses.I3Geometry()
  antGeometry.antennageo = dataclasses.I3AntennaGeoMap(GetGen1Antennas())

  print("Adding stations...")
  for ipos, pos in enumerate(stringLocations):
    AddStation(ipos+1, pos, scintGeometry.scintgeo, antGeometry.antennageo)

  print("Added {} panels".format(len(scintGeometry.scintgeo)))
  print("Added {} antennas".format(len(antGeometry.antennageo)))

  print("Adding DOMs...")
  #Only keep in-ice OMs
  omGeometry = frame['I3Geometry']
  newI3Geometry = dataclasses.I3Geometry()
  newI3Geometry.omgeo = dataclasses.I3OMGeoMap()
  for key in omGeometry.omgeo.keys():
    if key.om in [61, 62, 63, 64] and key.string <= 81: #Remove IT...
      if not args.keepit:
        continue

    newI3Geometry.omgeo[key] = omGeometry.omgeo[key]

  #Add the scint information
  outframe = icetray.I3Frame(icetray.I3Frame.Geometry) 

  outframe['I3Geometry'] = newI3Geometry
  outframe['I3ScintGeometry'] = scintGeometry
  outframe[radcube.GetDefaultGeometryName()] = antGeometry
  return outframe


##################################################
##################################################
##################################################

i3file = dataio.I3File(args.output, "w")

for frame in dataio.I3File(args.input,'r'):
  if frame.Stop == icetray.I3Frame.Geometry:

    outframe=ParseGeometryFrame(frame)
    i3file.push(outframe)
  else:
    i3file.push(frame)

print("Writing file", args.output)
i3file.close()
