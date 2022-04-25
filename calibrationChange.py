#!/usr/bin/env python
import numpy as np
from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
import os

gcd_fileBase = "/home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz"
gcd_fileSnow = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"


def ChangeGeometry(frame):
  #created newer version of geometry frame
  outframe  = icetray.I3Frame(icetray.I3Frame.Geometry)
  scintGeometry   = frame['I3ScintGeometry']
  antennaGeometry = frame['I3AntennaGeometry']
  ITGeometry      = frame['I3Geometry']  
  stageoOLD       = ITGeometry.stationgeo

  for f in dataio.I3File(gcd_fileSnow,'r'):
    if f.Stop == icetray.I3Frame.Geometry:
      geom = f['I3Geometry']
      stageo = geom.stationgeo
      for e,st in stageo:
        stageoOLD[e][0].snowheight = st[0].snowheight
        stageoOLD[e][1].snowheight = st[1].snowheight
        ##print(st[0].omkey_list) 
        ##print(stageoOLD[e], e)
        ##print(stageoOLD.keys())
  del frame['I3Geometry']
   
  outframe['I3ScintGeometry'] = scintGeometry
  outframe['I3AntennaGeometry'] = antennaGeometry
  outframe['I3Geometry'] = ITGeometry
  return outframe

# i3file = dataio.I3File("../simFiles/GCD-Survey-AntITScint_2020.02.24_Snow210305.i3.gz", "w")

# for frame in dataio.I3File(gcd_fileBase,'r'):
#   if frame.Stop == icetray.I3Frame.Geometry:
#     outframe=ChangeGeometry(frame)
#     i3file.push(outframe)
#   else:
#     i3file.push(frame)
# i3file.close()

def ChangeCalibration(frame):
  outframe = icetray.I3Frame(icetray.I3Frame.Calibration)
  GCDFile = dataio.I3File(gcd_fileSnow)
  cframe = GCDFile.pop_frame()
  while GCDFile.more() and not "I3Calibration" in cframe:
    cframe = GCDFile.pop_frame()
  return cframe


i3file = dataio.I3File("../simFiles/modified_GCD/GCD-Survey-AntITScint_2020.02.24_newCalibration.i3.gz", "w")

for frame in dataio.I3File(gcd_fileBase,'r'):
  if frame.Stop == icetray.I3Frame.Calibration:
    outframe=ChangeCalibration(frame)
    i3file.push(outframe)
  else:
    i3file.push(frame)
i3file.close()



