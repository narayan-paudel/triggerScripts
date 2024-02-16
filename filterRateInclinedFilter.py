#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio, filter_tools
from icecube.dataclasses import I3Constants
from icecube import topeventcleaning, tpx
from icecube import phys_services

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
mpl.use("Agg")
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mtick
from matplotlib.ticker import AutoMinorLocator,MultipleLocator
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap





plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)


import scipy.interpolate
import os



import glob

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--vertEvents", dest="vertEvents",
                    default=False, action="store_true", required=False,
                    help="plot vertical events")
args = parser.parse_args()

# usage python hitsInclinedShower.py --vertEvents
#to be run on envCurrent branch process2023


GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"

f = dataio.I3File(GCD)
dframe = f.pop_frame()
while f.more() and not "I3DetectorStatus" in dframe:
  dframe = f.pop_frame()
domStatus = dframe["I3DetectorStatus"].dom_status
om = dframe["I3Geometry"].omgeo
hg_doms = [ikey for ikey in om.keys() if (ikey.om in [61, 62, 63, 64] and ikey.string != 39 and domStatus[ikey].dom_gain_type == dataclasses.I3DOMStatus.High) or 
(ikey.om in [62, 63, 64] and ikey.string == 39 and domStatus[ikey].dom_gain_type == dataclasses.I3DOMStatus.High)]

# vertEvents = False
vertEvents = args.vertEvents
if not vertEvents:
  fileDir = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/"
  fileName = "combinedDeltaTIncl"
  plotSuffix = "Incl"
else:
  fileDir = "/home/enpaudel/dataExp/dataSetClean_VerticalLE/"
  fileName = "combinedDeltaTVert"
  plotSuffix = "Vert"

triggerListSelectDict = {"HG7_3":"IceTop7HG","HLC6":"IceTopSMT"}

fileDir = "/home/enpaudel/dataExp/dataSetCleanFilter/"
# fileList = sorted(glob.glob(fileDir+"*.i3.*"))[:10]
fileList = sorted(glob.glob(fileDir+"Fe*Evts.i3.*"))[:1]
fileName = "afterFilter.i3.gz"

print("fileList")

inclinationCut = 60 #degree
energyCut = 10**16 #eV
fileDir = "/home/enpaudel/dataExp/dataSetClean/"

exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]

outputDir = "/home/enpaudel/dataExp/dataSetCleanFilter/"
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

def triggerEfficiency(n_trig,n_total):
  # print("n_trig,n_total",n_trig,n_total)
  if n_total != 0:
    return (n_trig/n_total)
  else:
    return 0

def openingAngle(theta1,phi1,theta2,phi2):
  return np.arccos(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2))

def test7HG(frame):
  return frame["HG7_3"]>0


def excludeITSMT(frame):
  return frame["HLC6"]<1

def AddTotalCharge(frame):
  '''calculates total SLC charge in SLC tank pulses'''
  # print("Anything")
  if (frame.Has('OfflineIceTopSLCTankPulses')):
    # print("frame has OfflineIceTopSLCTankPulses")
    slc = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'OfflineIceTopSLCTankPulses')
    ITTotalChargeSLC = sum([sum([c.charge for c in slc[i]]) for i in slc.keys()])
    frame["ITTotalChargeSLC"] = dataclasses.I3Double(ITTotalChargeSLC)

def AddTotalCharge(frame,keys):
  '''calculates total SLC or HLC charges in tank pulses
  keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
  '''
  # print("Anything")
  for key in keys:    
    if (frame.Has(str(key))):
      # print("frame has", str(key))
      lc = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,str(key))
      ITTotalCharge = sum([sum([c.charge for c in lc[i]]) for i in lc.keys()])
      frame[str(key)+"TotalCharge"] = dataclasses.I3Double(ITTotalCharge)



def GetTimeToShowerPlane(particle, pos, time):
  t_plane = particle.GetTime() + (pos - particle.GetPos()).Dot(I3Position(particle.GetDir())) / I3Constants.c;
  return t_plane - time


def getSPTime(particle,omkey,pulseTime):
  '''
  returns a shower plane time.
  '''
  domgeo = geoFrame["I3Geometry"].omgeo[omkey]
  # pos_sc = radcube.GetShowerFromIC(domgeo.position - particle.pos, particle.dir)
  # t_offset = pos_sc.z/dataclasses.I3Constants.c
  # print("deltaT",pulseTime + t_offset)
  # return pulseTime + t_offset
  return pulseTime

def distanceToPlane(x,y,z,nx,ny,nz):
  '''
  when the plane pass through the origin
  point(x,y,z) to plane (nx,ny,nz)
  '''
  return (x*nx+y*ny+z*nz)/np.sqrt(nx**2+ny**2+ny**2)

def removePFrame(frame):
  return frame.Stop != icetray.I3Frame.Physics

def WilsonMean(nPass, nFail, z=1.0):  
  if nPass < 0 or nFail < 0:
      return 0
  return (nPass + z ** 2 / 2.0) / float(nPass + nFail + z ** 2)


def WilsonError(nPass, nFail, z=1.0):
  if nPass < 0 or nFail < 0:
      return 0

  if nPass + nFail == 0:
      return 0.5

  term = z / float(nPass + nFail + z ** 2)

  return term * (nPass * nFail / float(nPass + nFail) + z ** 2 / 4.0) ** 0.5

class showerObj(object):
  """docstring for showerObj"""
  def __init__(self):
    super(showerObj, self).__init__()
    # self.arg = arg
    



class TriggerEff(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.showerObjList = []


  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    ishw = showerObj()
    ishw.H4aWeight = frame["H4aWeight"].value
    ishw.incl_filt = frame["IceTopIncl_24"]
    ishw.HG7_3 = frame["HG7_3"]
    ishw.xcore = frame["MCPrimary"].pos.x
    ishw.ycore = frame["MCPrimary"].pos.y
    # print("cores",(xcore,ycore),(xcore_reco,ycore_reco),(xcore-xcore_reco,ycore-ycore_reco),r_diff)
    ishw.zenith_true = frame["MCPrimary"].dir.zenith
    ishw.azimuth_true = frame["MCPrimary"].dir.azimuth
    ishw.zenith_reco = frame["ShowerPlaneIncl"].dir.zenith
    ishw.azimuth_reco = frame["ShowerPlaneIncl"].dir.azimuth
    ishw.energy = frame["MCPrimary"].energy
    self.showerObjList.append(ishw)

  def plotTrigEfficiency(self,triggerType,containment):
    '''
    plots trigger efficiency in different zenith bins
    '''
    print("plotting trigger efficiency for ",triggerType)
    print("shower list",len(self.showerObjList))
    if containment == True:
      # evtList = containedEvents(evtList,640)
      evtList = [ishw for ishw in self.showerObjList if (ishw.xcore**2+ishw.ycore**2)<= 410**2]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    colorIter = iter(colorsCustom+colorsCustom)
    # colorIter = iter(colorsList)
    # energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBins = 10**np.linspace(5.0, 7.8, 29)
    for nbin, binStart in enumerate(sin2ZenBins[:-1]):
      lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
      highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
      evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith_true < highEdge]
      energyList = []
      efficiencyList = []
      wilsonMeanList = []
      errorListLow = []
      errorListHigh = []
      ncolor = colorsCustom2[nbin]
      for ebin, ebinStart in enumerate(energyBins[:-1]):
        lowEdge_E = energyBins[ebin]
        highEdge_E = energyBins[ebin+1]
        evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
        # totalEvts = len(evtEBin)
        weights = [ievt.H4aWeight for ievt in evtEBin]
        totalEvts = sum(weights)
        triggerList = [ievt.H4aWeight for ievt in evtEBin if abs(int(getattr(ievt,str(triggerType)).value)-1)<0.01]
        trigEvts = sum(triggerList)
        trigEff = triggerEfficiency(trigEvts,totalEvts)
        nPass = trigEff * len(weights)
        nFail = (1-trigEff) * len(weights)
        efficiencyList.append(trigEff)
        #####################################################################################
        #binomial interval
        # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68)
        #####################################################################################
        wilsonM = WilsonMean(nPass=nPass, nFail=nFail)
        wilsonErr = WilsonError(nPass=nPass, nFail=nFail)
        ErrorH = wilsonM + wilsonErr - trigEff
        ErrorL = trigEff - (wilsonM - wilsonErr)
        # print("wilson calculation",totalEvts,len(triggerList),trigEff,wilsonErr)
        wilsonMeanList.append(wilsonM)
        errorListLow.append(ErrorL)
        errorListHigh.append(ErrorH)
        energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))
      # ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt=markers2[ntrig],ls="-",lw=2.5,c=ncolor,label=str(triggerListSelectDict[itrigger]),alpha=1)
      ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"log$_{10}$ (E [eV])", fontsize=22)
    ax.set_ylabel(r"trigger efficiency", fontsize=22)
    # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
    ax.text(0.78,0.1,s=r"trig: {0}".format(triggerListSelectDict[triggerType]),size=13,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
    ax.text(0.78,0.2,s=r"snow: {0}".format("2021/03"),size=13,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
    # ax.set_xscale('log')
    ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
    ax.set_ylim(0,1.05)
    # ax.set_xlim(14,17)
    ax.set_xlim(14,16.8)
    # ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.grid(True,alpha=0.5)
    # l1=ax.legend(loc="upper left",fontsize=12)
    l1=ax.legend(loc="upper left",fontsize=12)
    point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"0.98")
    # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
    # ax.add_artist(l1)
    # ax.add_artist(l2)
    plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"EfficiencyFilt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"EfficiencyFilt.png",transparent=False,bbox_inches='tight')
    plt.close()

  def plotFilterEfficiency(self,filterType,containment):
    '''
    plots trigger efficiency in different zenith bins
    '''
    print("plotting trigger efficiency for ",triggerType)
    print("shower list",len(self.showerObjList))
    if containment == True:
      # evtList = containedEvents(evtList,640)
      evtList = [ishw for ishw in self.showerObjList if (ishw.xcore**2+ishw.ycore**2)<= 410**2]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    colorIter = iter(colorsCustom+colorsCustom)
    # colorIter = iter(colorsList)
    # energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBins = 10**np.linspace(5.0, 7.8, 29)
    for nbin, binStart in enumerate(sin2ZenBins[:-1]):
      lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
      highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
      evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith_true < highEdge]
      energyList = []
      efficiencyList = []
      wilsonMeanList = []
      errorListLow = []
      errorListHigh = []
      ncolor = colorsCustom2[nbin]
      for ebin, ebinStart in enumerate(energyBins[:-1]):
        lowEdge_E = energyBins[ebin]
        highEdge_E = energyBins[ebin+1]
        evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
        # totalEvts = len(evtEBin)
        weights = [ievt.H4aWeight for ievt in evtEBin]
        totalEvts = sum(weights)
        triggerList = [ievt.H4aWeight for ievt in evtEBin if abs(int(getattr(ievt,str(triggerType)).value)-1)<0.01]
        trigEvts = sum(triggerList)
        trigEff = triggerEfficiency(trigEvts,totalEvts)
        nPass = trigEff * len(weights)
        nFail = (1-trigEff) * len(weights)
        efficiencyList.append(trigEff)
        #####################################################################################
        #binomial interval
        # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68)
        #####################################################################################
        wilsonM = WilsonMean(nPass=nPass, nFail=nFail)
        wilsonErr = WilsonError(nPass=nPass, nFail=nFail)
        ErrorH = wilsonM + wilsonErr - trigEff
        ErrorL = trigEff - (wilsonM - wilsonErr)
        # print("wilson calculation",totalEvts,len(triggerList),trigEff,wilsonErr)
        wilsonMeanList.append(wilsonM)
        errorListLow.append(ErrorL)
        errorListHigh.append(ErrorH)
        energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))
      # ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt=markers2[ntrig],ls="-",lw=2.5,c=ncolor,label=str(triggerListSelectDict[itrigger]),alpha=1)
      ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"log$_{10}$ (E [eV])", fontsize=22)
    ax.set_ylabel(r"trigger efficiency", fontsize=22)
    # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
    ax.text(0.78,0.1,s=r"trig: {0}".format(triggerListSelectDict[triggerType]),size=13,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
    ax.text(0.78,0.2,s=r"snow: {0}".format("2021/03"),size=13,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
    # ax.set_xscale('log')
    ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
    ax.set_ylim(0,1.05)
    # ax.set_xlim(14,17)
    ax.set_xlim(14,16.8)
    # ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.grid(True,alpha=0.5)
    # l1=ax.legend(loc="upper left",fontsize=12)
    l1=ax.legend(loc="upper left",fontsize=12)
    point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"0.98")
    # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
    # ax.add_artist(l1)
    # ax.add_artist(l2)
    plt.savefig(plotFolder+"/trig"+str(filterType)+"cont"+str(containment)+"EfficiencyFilt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/trig"+str(filterType)+"cont"+str(containment)+"EfficiencyFilt.png",transparent=False,bbox_inches='tight')
    plt.close()

  def Finish(self):
    self.plotTrigEfficiency("HG7_3",containment=True)

    #####################################################
    #####################################################

name = ""

tray = I3Tray()
tray.AddModule("I3Reader","reader",
             # filenameList = args.input,
             filenameList = [GCD]+fileList,
             # filename = inFile,
            )

tray.AddModule(TriggerEff,"TEFF",
              # streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
              )


# tray.AddModule("I3Writer","i3writer",
#             filename=str(outputDir)+fileName,
#             Streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )

tray.Execute()
tray.Finish()