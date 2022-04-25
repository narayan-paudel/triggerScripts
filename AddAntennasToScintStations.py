#!/usr/bin/env python

# Reads in an existing GCD file with the optical strings and scintillator panels
# Antennas are added to the stations based on the scintillator locations
# It is expected that the panels at the center of the station have IDs 0 and 1

# ./AddAntennasToScintStations.py MyI3File.i3.gz --keepit

from icecube import icetray, dataclasses, dataio, radcube
from icecube.icetray import I3Units
import numpy as np
import os
from os.path import expandvars

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="Location of input GCD file to add stations to")
parser.add_argument("--output", type=str, default=ABS_PATH_HERE + "/GCD-AntennasAdded.i3.gz", help="Location of output GCD file after adding stations")
parser.add_argument("--keepit", action="store_true", help="If selected, IceTop will be kept in the file")
parser.add_argument("--scintperstn", type=int, default=8, help="Number of panels per station")
parser.add_argument("--antarmlength", type=float, default=72 / 2.0, help="Number of panels per station")
args = parser.parse_args()


def GetScintsFromGCD(filename):
    for frame in dataio.I3File(filename, "r"):
        if frame.Stop == icetray.I3Frame.Geometry:
            return frame["I3ScintGeometry"].scintgeo


def MakeAntennaGeo(name, position, heightAboveSnow=1.2):
    """
  Makes a single instance of a antenna geometry item
  Arguments:
  name - the name (str) of the antenna (don't really matter)
  position - I3Position of the antenna location
  heightAboveSnow - how high above the snow the antennas is

  Returns:
  I3AntennaGeo instance corresponding to the inputs
  """
    antennaGeo = dataclasses.I3AntennaGeo()
    antennaGeo.position = position
    antennaGeo.orientation = dataclasses.I3Orientation(1, 0, 0, 0, 0, 1)
    antennaGeo.heightAboveSnow = heightAboveSnow
    antennaGeo.antennaName = name
    antennaGeo.antennatype = dataclasses.I3AntennaGeo.AntennaType.SKALA2
    antennaGeo.cableLength = 50 * I3Units.m

    return antennaGeo


def AddAntennasToStation(stationId, scintgeomap, antennageomap):
    """
    Given a station ID corresponding to a station in the scintgeomap,
    the panels in the scint map are analyzed and an antenna is added to
    the antennageomap for each spoke. The antenna is placed a distance
    down the spoke corresponding to the distance given in args.

    Arguments:
    stationId - int corresponding to the station number (used for keys)
    scintgeomap - I3ScintGeoMap which is filled with scint events
    antennageomap - I3AntennaGeoMap to be added to by this function
    """

    scintkeys = [key for key in scintgeomap.keys() if key.station == stationId]
    assert len(scintkeys) == args.scintperstn

    stationCenter = scintgeomap[dataclasses.ScintKey(stationId, 0)].position
    stationCenter += scintgeomap[dataclasses.ScintKey(stationId, 1)].position
    stationCenter /= 2.0

    antPos = []

    nSpokes = int((args.scintperstn - 2) / 2)
    for ispoke in range(nSpokes):
        panel1 = scintgeomap[dataclasses.ScintKey(stationId, 2 + ispoke * 2)].position
        panel2 = scintgeomap[dataclasses.ScintKey(stationId, 3 + ispoke * 2)].position

        posSpoke = (panel1 + panel2) / 2.0

        vec = posSpoke - stationCenter
        vec.normalize()

        pos = vec * args.antarmlength + stationCenter
        antkey = dataclasses.AntennaKey(stationId, ispoke + 1)
        antennageomap[antkey] = MakeAntennaGeo("stn{}_ant{}".format(stationId, ispoke + 1), pos)


def ParseGeometryFrame(gFrame):
    """
    Reads in the info of the base GCD file, adds stations on top of each of the
    non-Gen1 strings. Adds all the stations from the Gen1 array. Adds all the
    DOMS back into the GFrame. Returns a new copy of the frame with updated info.

    Arguments:
    gFrame - Geometry frame from the baseline DOM layout

    Returns:
    Geometry frame with all the additional fields
    """

    assert frame.Stop == icetray.I3Frame.Geometry  # Make sure it is a G-Frame

    scintGeometry = dataclasses.I3Geometry()
    scintGeometry.scintgeo = dataclasses.I3ScintGeoMap(GetScintsFromGCD(args.input))

    scintStations = list(dict.fromkeys([key.station for key in scintGeometry.scintgeo.keys()]))
    scintStations.sort()

    antGeometry = dataclasses.I3Geometry()
    antGeometry.antennageo = dataclasses.I3AntennaGeoMap()

    print("Adding antennas stations...")
    for stationId in scintStations:
        AddAntennasToStation(stationId, scintGeometry.scintgeo, antGeometry.antennageo)

    print("Added {} panels".format(len(scintGeometry.scintgeo)))
    print("Added {} antennas".format(len(antGeometry.antennageo)))

    print("Adding DOMs...")
    # Only keep in-ice OMs
    omGeometry = frame["I3Geometry"]
    newI3Geometry = dataclasses.I3Geometry()
    newI3Geometry.omgeo = dataclasses.I3OMGeoMap()
    for key in omGeometry.omgeo.keys():
        if key.om in [61, 62, 63, 64] and key.string <= 81:  # Remove IT...
            if not args.keepit:
                continue

        newI3Geometry.omgeo[key] = omGeometry.omgeo[key]

    if args.keepit:
        newI3Geometry.stationgeo = omGeometry.stationgeo

    # Add the scint information
    outframe = icetray.I3Frame(icetray.I3Frame.Geometry)

    outframe["I3Geometry"] = newI3Geometry
    outframe["I3ScintGeometry"] = scintGeometry
    outframe[radcube.GetDefaultGeometryName()] = antGeometry
    return outframe


##################################################
##################################################
##################################################

i3file = dataio.I3File(args.output, "w")

for frame in dataio.I3File(args.input, "r"):
    if frame.Stop == icetray.I3Frame.Geometry:

        outframe = ParseGeometryFrame(frame)
        i3file.push(outframe)
    else:
        i3file.push(frame)

print("Writing file", args.output)
i3file.close()
