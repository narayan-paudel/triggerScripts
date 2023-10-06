#!/usr/bin/env python3

'''
This example shows how to include IceTop DOMSets to the GCD file in addition to default InIce DOMSets
'''
from icecube.icetray import I3Tray, I3Units
from icecube import dataclasses,icetray,dataio
from icecube.trigger_sim import GetDefaultDOMSets
from icecube.trigger_sim.InjectDefaultDOMSets import InjectDefaultDOMSets


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input',"-i",type=str,default="",help="input simulation GCD for IceTop")
parser.add_argument('--output',"-o",type=str,default="testGCD.i3.gz",help='output simulation GCD for IceTop')
args = parser.parse_args()

icetop_domset_definitions = {}

def IceTopDoms(omkey, omgeo):
  """selects IceTop doms (domset 3)"""
  return (omkey.string >= 1 and omkey.string <= 81 and omkey.om >= 61 and omkey.om <= 64 )

def highGainInFillIceTopDoms(omkey, omgeo):
  """selects high gain InFill IceTop doms (domset 8)"""
  DOMSET_8_IT_STATIONS = [26,36,46,79,80,81]
  return (omkey.string in DOMSET_8_IT_STATIONS and 
    (omkey.om == (61 + (omkey.string==26)) or
      omkey.om == 63))

def highGainIceTopDoms(omkey, omgeo):
  """selects high gain IceTop doms (domset 12)"""
  return (omkey.string >= 1 and omkey.string <= 81 and 
    (omkey.om == (61 + (omkey.string==26 or omkey.string==39 or omkey.string==74)) or
      (omkey.om == (63 + (omkey.string==67)))))

icetop_domset_definitions[3] = IceTopDoms
icetop_domset_definitions[8] = highGainInFillIceTopDoms
icetop_domset_definitions[12] = highGainIceTopDoms

tray = I3Tray()
tray.AddModule("I3Reader","reader",
             filename=args.input,
            )

tray.Add(InjectDefaultDOMSets,
         NewDefinitionMap = icetop_domset_definitions)

tray.AddModule("I3Writer","i3writer",
            filename=args.output,
            streams=[icetray.I3Frame.Geometry,icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus],
            )

tray.Execute()
tray.Finish()

