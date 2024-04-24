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
# fileList = sorted(glob.glob(fileDir+"*.i3.*"))
fileList = sorted(glob.glob(fileDir+"Fe*Evts.i3.*"))[:10]
fileName = "misReco.i3.gz"

print("fileList")

inclinationCut = 60 #degree
energyCut = 10**16 #eV
fileDir = "/home/enpaudel/dataExp/dataSetClean/"

exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]

outputDir = "/home/enpaudel/dataExp/dataSetCleanFilterMisReco/"
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

def triggerEfficiency(n_trig,n_total):
  # print("n_trig,n_total",n_trig,n_total)
  if n_total != 0:
    return (n_trig/n_total)
  else:
    return 0

def filterEfficiency(n_trig,n_total):
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
    psm_HLC = frame["OfflineIceTopHLCTankPulsesCleanTimeCleanCharge"]
    if psm_HLC.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm_HLC.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm_HLC = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, "OfflineIceTopHLCTankPulsesCleanTimeCleanCharge")
    ishw.nHLC = len(psm_HLC)
    ishw.Qtot_HLC = sum([sum([c.charge for c in psm_HLC[i]]) for i in psm_HLC.keys()])
    psm_SLC = frame["OfflineIceTopSLCTankPulsesCleanTimeCleanCharge"]
    if psm_SLC.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm_SLC.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm_SLC = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, "OfflineIceTopSLCTankPulsesCleanTimeCleanCharge")
    ishw.nSLC = len(psm_SLC)
    ishw.Qtot_SLC = sum([sum([c.charge for c in psm_SLC[i]]) for i in psm_SLC.keys()])
    ishw.H4aWeight = frame["H4aWeight"].value
    ishw.IceTopIncl_24 = frame["IceTopIncl_24"]
    ishw.HG7_3 = frame["HG7_3"]
    ishw.HLC6 = frame["HLC6"]
    ishw.xcore = frame["MCPrimary"].pos.x
    ishw.ycore = frame["MCPrimary"].pos.y
    ishw.xcoreCOG = frame["ShowerCOGIncl"].pos.x
    ishw.ycoreCOG = frame["ShowerCOGIncl"].pos.y
    # print("cores",(xcore,ycore),(xcore_reco,ycore_reco),(xcore-xcore_reco,ycore-ycore_reco),r_diff)
    ishw.zenith_true = frame["MCPrimary"].dir.zenith
    # print("zenith",frame["MCPrimary"].dir.zenith,np.rad2deg(frame["MCPrimary"].dir.zenith))
    ishw.azimuth_true = frame["MCPrimary"].dir.azimuth
    ishw.zenith_reco = frame["ShowerPlaneIncl"].dir.zenith
    ishw.zenith_lapu = frame["Laputop"].dir.zenith
    ishw.azimuth_reco = frame["ShowerPlaneIncl"].dir.azimuth
    ishw.azimuth_lapu = frame["Laputop"].dir.azimuth
    ishw.energy = frame["MCPrimary"].energy
    frame["zenith_true"] = dataclasses.I3Double(np.rad2deg(frame["MCPrimary"].dir.zenith))
    frame["azi_true"] = dataclasses.I3Double(np.rad2deg(frame["MCPrimary"].dir.azimuth))
    frame["zenith_reco"] = dataclasses.I3Double(np.rad2deg(frame["ShowerPlaneIncl"].dir.zenith))
    frame["azi_reco"] = dataclasses.I3Double(np.rad2deg(frame["ShowerPlaneIncl"].dir.azimuth))
    self.PushFrame(frame)
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
    print("plotting filter efficiency for ",filterType)
    print("shower list",len(self.showerObjList))
    triggerType = "HG7_3"
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
        filterList = [ievt.H4aWeight for ievt in evtEBin if abs(int(getattr(ievt,str(filterType)).value)-1)<0.01 and abs(int(getattr(ievt,str(triggerType)).value)-1)<0.01]
        trigEvts = sum(triggerList)
        filtEvts = sum(filterList)
        # trigEff = triggerEfficiency(trigEvts,totalEvts)
        filtEff = filterEfficiency(filtEvts,trigEvts)
        nPass = filtEff * len(weights)
        nFail = (1-filtEff) * len(weights)
        efficiencyList.append(filtEff)
        #####################################################################################
        #binomial interval
        # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68)
        #####################################################################################
        wilsonM = WilsonMean(nPass=nPass, nFail=nFail)
        wilsonErr = WilsonError(nPass=nPass, nFail=nFail)
        ErrorH = wilsonM + wilsonErr - filtEff
        ErrorL = filtEff - (wilsonM - wilsonErr)
        # print("wilson calculation",totalEvts,len(triggerList),trigEff,wilsonErr)
        wilsonMeanList.append(wilsonM)
        errorListLow.append(ErrorL)
        errorListHigh.append(ErrorH)
        energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))
      # ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt=markers2[ntrig],ls="-",lw=2.5,c=ncolor,label=str(triggerListSelectDict[itrigger]),alpha=1)
      ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"log$_{10}$ (E [eV])", fontsize=22)
    ax.set_ylabel(r"filter efficiency", fontsize=22)
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
    plt.savefig(plotFolder+"/filt"+str(filterType)+"cont"+str(containment)+"EfficiencyFilt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/filt"+str(filterType)+"cont"+str(containment)+"EfficiencyFilt.png",transparent=False,bbox_inches='tight')
    plt.close()


  def plotTrigEfficiencyZenith(self,triggerType,containment,wilson,weighting):
    '''
    plots trigger efficiency in different zenith bins
    '''
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    print("plotting trigger efficiency for ",triggerType)
    if containment == True:
      # evtList = containedEvents(evtList,640)
      evtList = [ishw for ishw in self.showerObjList if (ishw.xcore**2+ishw.ycore**2)<= 410**2]
    # sin2ZenBins = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.822]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    colorIter = iter(colorsCustom+colorsCustom)
    # colorIter = iter(colorsList)
    # zorders = [0,1,2,3,4,5,6]
    zorders = [60,50,40,30,20,10,9]
    for nbin, binStart in enumerate(energyBinsShort[:-1]):
      lowEdge = energyBinsShort[nbin]
      highEdge = energyBinsShort[nbin+1]
      evtEnergyBin = [ievt for ievt in evtList if lowEdge <= ievt.energy < highEdge]
      zenithList = []
      efficiencyList = []
      wilsonMeanList = []
      errorListLow = []
      errorListHigh = []
      ncolor = colorsCustom2[nbin]
      for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
        lowEdge_Z = sin2ZenBins[ebin]
        highEdge_Z = sin2ZenBins[ebin+1]
        evtZBin = [ievt for ievt in evtEnergyBin if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
        weights = [ievt.H4aWeight for ievt in evtZBin]
        totalEvts = sum(weights)
        triggerList = [ievt.H4aWeight for ievt in evtZBin if abs(getattr(ievt,triggerType).value-1)<0.01]
        trigEvts = sum(triggerList)
        trigEff = triggerEfficiency(trigEvts,totalEvts)
        nPass = trigEff * len(weights)
        nFail = (1-trigEff) * len(weights)
        #########################################################
        # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68) #binomial interval
        ##########################################################
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
        efficiencyList.append(trigEff)
        zenithList.append((np.arcsin(np.sqrt(lowEdge_Z))+np.arcsin(np.sqrt(highEdge_Z)))/2.0*180/np.pi)
      # ax.plot(zenithList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
      # ax.errorbar(zenithList,wilsonMeanList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
      # ax.errorbar(zenithList,wilsonMeanList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="none",ls = 2.5,c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
      ax.errorbar(zenithList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",ls="-",lw = 2.5,c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
      # ax.plot(zenithList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,zorder=zorders[nbin],alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"zenith [$^{\circ}$]", fontsize=22)
    ax.set_ylabel(r"trigger efficiency", fontsize=22)
    ax.text(0.02,0.62,s="IceCube Preliminary",color="red",size=14,fontweight='bold',horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
    # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
    ax.text(0.02,0.73,s=r"trig: {0}".format(triggerListSelectDict[triggerType]),size=14,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,zorder=10002)
    ax.text(0.02,0.67,s=r"snow: {0}".format("2021/03"),size=14,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,zorder=10001)
    # ax.set_xscale('log')
    ax.set_xticks(np.linspace(0,70,8))
    yline=0.98
    ax.axhline(y=yline,xmin=0,xmax=1,color="gray",linestyle="--",zorder=1,lw=2.0)
    ax.set_ylim(0,1.05)
    ax.set_xlim(0,65)
    # ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.grid(True,alpha=0.5)
    l1=ax.legend(loc="lower left",title=r"log$_{10}$ [E/eV]",ncol=1,fontsize=18)
    l1.set_zorder(10000)
    # l1=ax.legend(loc="upper left",title=r"log$_{10}$ [E/eV]",title_fontsize=18,fontsize=18).set_zorder(10000)
    # l1=ax.legend(loc="upper left",fontsize=18).set_zorder(10000)
    # plt.setp(l1.get_title(),fontsize=18)
    # plt.legend.set_title(r"log$_{10}$ (E/eV)",prop={'size':14})
    # point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"{yline:.2f}")
    point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='gray', marker='',markersize=5, label=r"0.98")
    l2 = ax.legend(handles=[point_dash],loc=(0.02,0.76),fontsize=18,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
    ax.add_artist(l1)
    ax.add_artist(l2)
    plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"Weight"+str(weighting)+"EfficiencyZenithFilt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"Weight"+str(weighting)+"EfficiencyZenithFilt.png",transparent=False,bbox_inches='tight')
    plt.close()


  def plotFilterEfficiencyZenith(self,filterType,containment,wilson,weighting):
    '''
    plots trigger efficiency in different zenith bins
    '''
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    triggerType = "HG7_3"
    if containment == True:
      # evtList = containedEvents(evtList,640)
      evtList = [ishw for ishw in self.showerObjList if (ishw.xcore**2+ishw.ycore**2)<= 410**2]
    # sin2ZenBins = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.822]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    colorIter = iter(colorsCustom+colorsCustom)
    # colorIter = iter(colorsList)
    # zorders = [0,1,2,3,4,5,6]
    zorders = [60,50,40,30,20,10,9]
    for nbin, binStart in enumerate(energyBinsShort[:-1]):
      lowEdge = energyBinsShort[nbin]
      highEdge = energyBinsShort[nbin+1]
      evtEnergyBin = [ievt for ievt in evtList if lowEdge <= ievt.energy < highEdge]
      zenithList = []
      efficiencyList = []
      wilsonMeanList = []
      errorListLow = []
      errorListHigh = []
      ncolor = colorsCustom2[nbin]
      for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
        lowEdge_Z = sin2ZenBins[ebin]
        highEdge_Z = sin2ZenBins[ebin+1]
        evtZBin = [ievt for ievt in evtEnergyBin if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
        weights = [ievt.H4aWeight for ievt in evtZBin]
        totalEvts = sum(weights)
        triggerList = [ievt.H4aWeight for ievt in evtZBin if abs(getattr(ievt,triggerType).value-1)<0.01]
        filterList = [ievt.H4aWeight for ievt in evtZBin if abs(getattr(ievt,triggerType).value-1)<0.01 and abs(getattr(ievt,filterType).value-1)<0.01]
        trigEvts = sum(triggerList)
        filtEvts = sum(filterList)
        trigEff = triggerEfficiency(filtEvts,trigEvts)
        nPass = trigEff * len(weights)
        nFail = (1-trigEff) * len(weights)
        #########################################################
        # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68) #binomial interval
        ##########################################################
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
        efficiencyList.append(trigEff)
        zenithList.append((np.arcsin(np.sqrt(lowEdge_Z))+np.arcsin(np.sqrt(highEdge_Z)))/2.0*180/np.pi)
      # ax.plot(zenithList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
      # ax.errorbar(zenithList,wilsonMeanList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
      # ax.errorbar(zenithList,wilsonMeanList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="none",ls = 2.5,c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
      ax.errorbar(zenithList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",ls="-",lw = 2.5,c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
      # ax.plot(zenithList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,zorder=zorders[nbin],alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"zenith [$^{\circ}$]", fontsize=22)
    ax.set_ylabel(r"filter efficiency", fontsize=22)
    ax.text(0.02,0.62,s="IceCube Preliminary",color="red",size=14,fontweight='bold',horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
    # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
    ax.text(0.02,0.73,s=r"trig: {0}".format(triggerListSelectDict[triggerType]),size=14,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,zorder=10002)
    ax.text(0.02,0.67,s=r"snow: {0}".format("2021/03"),size=14,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,zorder=10001)
    # ax.set_xscale('log')
    ax.set_xticks(np.linspace(0,70,8))
    yline=0.98
    ax.axhline(y=yline,xmin=0,xmax=1,color="gray",linestyle="--",zorder=1,lw=2.0)
    ax.set_ylim(0,1.05)
    ax.set_xlim(0,65)
    # ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.grid(True,alpha=0.5)
    l1=ax.legend(loc="lower left",title=r"log$_{10}$ [E/eV]",ncol=1,fontsize=18)
    l1.set_zorder(10000)
    # l1=ax.legend(loc="upper left",title=r"log$_{10}$ [E/eV]",title_fontsize=18,fontsize=18).set_zorder(10000)
    # l1=ax.legend(loc="upper left",fontsize=18).set_zorder(10000)
    # plt.setp(l1.get_title(),fontsize=18)
    # plt.legend.set_title(r"log$_{10}$ (E/eV)",prop={'size':14})
    # point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"{yline:.2f}")
    point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='gray', marker='',markersize=5, label=r"0.98")
    l2 = ax.legend(handles=[point_dash],loc=(0.02,0.76),fontsize=18,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
    ax.add_artist(l1)
    ax.add_artist(l2)
    plt.savefig(plotFolder+"/filt"+str(filterType)+"cont"+str(containment)+"Weight"+str(weighting)+"EfficiencyZenith.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/filt"+str(filterType)+"cont"+str(containment)+"Weight"+str(weighting)+"EfficiencyZenith.png",transparent=False,bbox_inches='tight')
    plt.close()


  def plotEnergyFlux(self,triggerType,yscale,suffix,energyScale,containment):
    """
    plots energy flux
    """
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    eventList = self.showerObjList
    if containment == True:
      # evtList = containedEvents(evtList,640)
      eventList = [ishw for ishw in eventList if (ishw.xcore**2+ishw.ycore**2)<= 410**2]
      # eventList = containedEvents(eventList,800)
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType).value-0)>0.01]
    # eventList = [ievt for ievt in eventList if 10**6.9 <= ievt.energy < 10**8.0]    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    colorIter = iter(colorsCustom+colorsCustom)
    hitBins = np.linspace(14.0,17.0,31)
    # hitBins = np.linspace(14.0,16.8,29)
    totalRate = 0
    for nbin, binStart in enumerate(sin2ZenBins[:-1]):
      lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
      highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
      evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith_true < highEdge]
      energyList = []
      neventList = []
      # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
      ncolor = next(colorIter)
      ax,histSum = self.plotSteps(evtZenBin,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),energyScale,ncolor=colorsCustom2[nbin])
      totalRate += histSum
    # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
    # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"log$_{10}$ (E [eV])", fontsize=22)
    # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
    ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
    ax.set_yscale(yscale) 
    ax.set_ylim(10**-5,10**0.0)
    ax.set_xlim(14.0,17.0)
    ax.set_xlim(14.0,16.8)
    # ax.set_xscale('log')
    ax.grid(True,alpha=0.7)
    ax.text(0.82,0.82,s=r"{0} rate:{1:.1f} Hz".format(triggerListSelectDict[triggerType],totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
    # ax.set_title(key,fontsize=16)
    ax.legend(fontsize=10,ncol=3,loc="lower center")
    plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+"Filt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+"Filt.png",transparent=False,bbox_inches='tight')
    plt.close()
  def plotSteps(self,triggeredEvts,ax,legendLabel,energyScale,ncolor):
    hitBins = np.linspace(14.0,17.0,31)
    energy = [ievt.energy for ievt in triggeredEvts]
    weights = [ievt.H4aWeight for ievt in triggeredEvts]
    # weights_direct = [ievt.directWeight for ievt in triggeredEvts]
    energy = np.log10(energy)+9.0
    hist,binEdge = np.histogram(energy,hitBins,weights=[w for w in weights])
    # binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
    # print("sum of simulated rate",sum(hist))
    # if "total" in legendLabel:
      # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
    # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
      # ax.set_ylim(10**-4,20)
    binCenter = (binEdge[:-1]+binEdge[1:])/2.0
    if str(energyScale) == "0.0":
      H = [h for h in hist]   
    else:
      H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
    # ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
    # ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.6f} Hz".format(sum(hist)),color=ncolor,alpha=1)
    ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.3f} Hz".format(sum(hist)),color=ncolor,alpha=1)
    return ax,sum(hist)

  def plotZenithFlux(self,triggerType,yscale,suffix,containment):
    """
    plots energy flux
    """
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    eventList = self.showerObjList
    if containment == True:
      # evtList = containedEvents(evtList,640)
      eventList = [ishw for ishw in eventList if (ishw.xcore**2+ishw.ycore**2)<= 410**2]
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType).value-1)<0.01]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    colorIter = iter(colorsCustom+colorsCustom)
    hitBins = np.linspace(14.0,17.0,31)
    totalRate = 0
    for nbin, binStart in enumerate(energyBinsShort[:-1]):
      lowEdge = energyBinsShort[nbin]
      highEdge = energyBinsShort[nbin+1]
      evtEnergyBin = [ievt for ievt in eventList if lowEdge <= ievt.energy < highEdge]
      zenithList = []
      neventList = []
      # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
      ncolor = next(colorIter)
      ax,histSum = self.plotZenithSteps(evtEnergyBin,ax,r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),ncolor=colorsCustom2[nbin])
      totalRate += histSum
    # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
    # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"zenith [$^{\circ}$]", fontsize=22)
    ax.set_xticks(np.linspace(0,70,8))
    # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
    ax.set_ylabel(r"log10(rate [Hz])", fontsize=20)
    # ax.set_yscale(yscale) 
    # ax.set_ylim(0,65)
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.set_xlim(0,65)
    ax.set_ylim(-4.5,0.5)
    # ax.set_xscale('log')
    ax.grid(True,alpha=0.7)
    ax.text(0.82,0.82,s=r"{0} rate:{1:.2f} Hz".format(triggerListSelectDict[triggerType],totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
    # ax.set_title(key,fontsize=16)
    # ax.legend(fontsize=10,ncol=3,loc="lower center",columnspacing=0.4)
    ax.legend(fontsize=10,ncol=1,loc="upper left",columnspacing=0.4)
    plt.savefig(plotFolder+"/ZenithSpecTrig"+str(triggerListSelectDict[triggerType])+str(suffix)+"scaleCont"+str(containment)+"Filt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/ZenithSpecTrig"+str(triggerListSelectDict[triggerType])+str(suffix)+"scaleCont"+str(containment)+"Filt.png",transparent=False,bbox_inches='tight')
    plt.close()

  def plotZenithSteps(self,triggeredEvts,ax,legendLabel,ncolor):
    zenithBins = [np.arcsin(np.sqrt(izen)) for izen in sin2ZenBins]
    zenith = [ievt.zenith_true for ievt in triggeredEvts]
    weights = [ievt.H4aWeight for ievt in triggeredEvts]
    # weights_direct = [ievt.directWeight for ievt in triggeredEvts]
    hist,binEdge = np.histogram(zenith,zenithBins,weights=[w for w in weights])
    # binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
    # print("sum of simulated rate",sum(hist))
    # if "total" in legendLabel:
      # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
    # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
      # ax.set_ylim(10**-4,20)
    # print("binedge",(binEdge[:-1]+binEdge[1:])/2.0)
    binCenter = (binEdge[:-1]+binEdge[1:])/2.0*180.0/np.pi
    H = [np.log10(h) for h in hist]   
    # ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
    ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.4f} Hz".format(sum(hist)),color=ncolor,alpha=1)
    return ax,sum(hist)


  def plotEnergyFluxFilter(self,triggerType,filterType,yscale,suffix,energyScale,containment):
    """
    plots energy flux
    """
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    eventList = self.showerObjList
    if containment == True:
      # evtList = containedEvents(evtList,640)
      eventList = [ishw for ishw in eventList if (ishw.xcore**2+ishw.ycore**2)<= 410**2]
      # eventList = containedEvents(eventList,800)
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType).value-0)>0.01 and abs(getattr(ievt,filterType).value-0)>0.01]
    # eventList = [ievt for ievt in eventList if 10**6.9 <= ievt.energy < 10**8.0]    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    colorIter = iter(colorsCustom+colorsCustom)
    hitBins = np.linspace(14.0,17.0,31)
    # hitBins = np.linspace(14.0,16.8,29)
    totalRate = 0
    for nbin, binStart in enumerate(sin2ZenBins[:-1]):
      lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
      highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
      evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith_true < highEdge]
      energyList = []
      neventList = []
      # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
      ncolor = next(colorIter)
      ax,histSum = self.plotSteps(evtZenBin,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),energyScale,ncolor=colorsCustom2[nbin])
      totalRate += histSum
    # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
    # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"log$_{10}$ (E [eV])", fontsize=22)
    # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
    ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
    ax.set_yscale(yscale) 
    ax.set_ylim(10**-5,10**0.0)
    ax.set_xlim(14.0,17.0)
    ax.set_xlim(14.0,16.8)
    # ax.set_xscale('log')
    ax.grid(True,alpha=0.7)
    ax.text(0.82,0.82,s=r"{0} rate:{1:.1f} Hz".format(triggerListSelectDict[triggerType],totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
    # ax.set_title(key,fontsize=16)
    ax.legend(fontsize=10,ncol=3,loc="lower center")
    plt.savefig(plotFolder+"/energySpecFilt"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+"Filt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/energySpecFilt"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+"Filt.png",transparent=False,bbox_inches='tight')
    plt.close()

  def plotZenithFluxFilter(self,triggerType,filterType,yscale,suffix,containment):
    """
    plots energy flux
    """
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    eventList = self.showerObjList
    if containment == True:
      # evtList = containedEvents(evtList,640)
      eventList = [ishw for ishw in eventList if (ishw.xcore**2+ishw.ycore**2)<= 410**2]
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType).value-1)<0.01 and abs(getattr(ievt,filterType).value-1)<0.01]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(ncols=1,nrows=1)
    ax = fig.add_subplot(gs[0])
    colorIter = iter(colorsCustom+colorsCustom)
    hitBins = np.linspace(14.0,17.0,31)
    totalRate = 0
    for nbin, binStart in enumerate(energyBinsShort[:-1]):
      lowEdge = energyBinsShort[nbin]
      highEdge = energyBinsShort[nbin+1]
      evtEnergyBin = [ievt for ievt in eventList if lowEdge <= ievt.energy < highEdge]
      zenithList = []
      neventList = []
      # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
      ncolor = next(colorIter)
      ax,histSum = self.plotZenithSteps(evtEnergyBin,ax,r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),ncolor=colorsCustom2[nbin])
      totalRate += histSum
    # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
    # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.set_xlabel(r"zenith [$^{\circ}$]", fontsize=22)
    ax.set_xticks(np.linspace(0,70,8))
    # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
    ax.set_ylabel(r"log10(rate [Hz])", fontsize=20)
    # ax.set_yscale(yscale) 
    # ax.set_ylim(0,65)
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.set_xlim(0,65)
    ax.set_ylim(-4.5,0.5)
    # ax.set_xscale('log')
    ax.grid(True,alpha=0.7)
    ax.text(0.82,0.82,s=r"{0} rate:{1:.2f} Hz".format(triggerListSelectDict[triggerType],totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
    # ax.set_title(key,fontsize=16)
    # ax.legend(fontsize=10,ncol=3,loc="lower center",columnspacing=0.4)
    ax.legend(fontsize=10,ncol=1,loc="upper left",columnspacing=0.4)
    plt.savefig(plotFolder+"/ZenithSpecFilt"+str(triggerListSelectDict[triggerType])+str(suffix)+"scaleCont"+str(containment)+"Filt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/ZenithSpecFilt"+str(triggerListSelectDict[triggerType])+str(suffix)+"scaleCont"+str(containment)+"Filt.png",transparent=False,bbox_inches='tight')
    plt.close()
  
  def plotHist(self,triggerType,filterType):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(0,1,21)
    zenithAll_lowE = [np.sin(ishwr.zenith_true)**2 for ishwr in self.showerObjList if ishwr.energy*10**9<10**16]
    zenithAll_highE = [np.sin(ishwr.zenith_true)**2 for ishwr in self.showerObjList if ishwr.energy*10**9>=10**16]
    zenithTrigger_lowE = [np.sin(ishwr.zenith_true)**2 for ishwr in self.showerObjList if abs(getattr(ishwr,triggerType).value-1)<0.01 and ishwr.energy*10**9<10**16]
    zenithTrigger_highE = [np.sin(ishwr.zenith_true)**2 for ishwr in self.showerObjList if abs(getattr(ishwr,triggerType).value-1)<0.01 and ishwr.energy*10**9>=10**16]
    zenithFilter_lowE = [np.sin(ishwr.zenith_true)**2 for ishwr in self.showerObjList if abs(getattr(ishwr,triggerType).value-1)<0.01 and abs(getattr(ishwr,filterType).value-1)<0.01 and ishwr.energy*10**9<10**16]
    zenithFilter_highE = [np.sin(ishwr.zenith_true)**2 for ishwr in self.showerObjList if abs(getattr(ishwr,triggerType).value-1)<0.01 and abs(getattr(ishwr,filterType).value-1)<0.01 and ishwr.energy*10**9>=10**16]
    ax.hist(zenithAll_lowE,bins=bins,histtype="step",label="All evts,E<10PeV")
    ax.hist(zenithAll_highE,bins=bins,histtype="step",label="All evts,E>=10PeV")
    ax.hist(zenithTrigger_lowE,bins = bins,histtype="step",label="HG7_3 trig,E<10PeV")
    ax.hist(zenithTrigger_highE,bins = bins,histtype="step",label="HG7_3 trig,E>=10PeV")
    ax.hist(zenithFilter_lowE,bins=bins,histtype="step",label="IceTopIncl_24 filt,E<10PeV")
    ax.hist(zenithFilter_highE,bins=bins,histtype="step",label="IceTopIncl_24 filt,E>=10PeV")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=11)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("count", fontsize=18)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel(r"sin$^2$$\theta$", fontsize=18)
    # ax.set_yscale("log")
    ax.legend(fontsize="11")
    plt.savefig(plotFolder+"/hist"+str(filterType)+".pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/hist"+str(filterType)+".png",transparent=False,bbox_inches='tight')

  def plotHistWeighted(self,triggerType,filterType):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(0,1,21)
    zenithAll = [np.sin(ishwr.zenith_true)**2 for ishwr in self.showerObjList]
    zenithAll_reco = [np.sin(ishwr.zenith_reco)**2 for ishwr in self.showerObjList]
    zenithAll_lapu = [np.sin(ishwr.zenith_lapu)**2 for ishwr in self.showerObjList]
    weightAll = [ishwr.H4aWeight for ishwr in self.showerObjList]
    ax.hist(zenithAll,bins=bins,weights=weightAll,histtype="step",label="MC true")
    ax.hist(zenithAll_reco,bins=bins,weights=weightAll,histtype="step",label="plane")
    ax.hist(zenithAll_lapu,bins=bins,weights=weightAll,histtype="step",label="Laputop")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=11)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("count", fontsize=18)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel(r"sin$^2$$\theta$", fontsize=18)
    # ax.set_yscale("log")
    ax.legend(fontsize="11")
    plt.savefig(plotFolder+"/hist"+str(filterType)+"Weight.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/hist"+str(filterType)+"Weight.png",transparent=False,bbox_inches='tight')

  def plotHistWeightedTest(self,triggerType,filterType,containment,showerObjList):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    zenCut = 45 #in degrees    
    sin_zenith = False
    if containment == True:
      showerObjList = [ishwr for ishwr in showerObjList if (ishwr.xcore**2+ishwr.ycore**2)<= 410**2]
    # showerObjList = [ishwr for ishwr in self.showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45 and np.rad2deg(ishwr.zenith_lapu) < 45]
    showerObjList = [ishwr for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= zenCut]
    # for ishwr in showerObjList:
    #   print("zenith",np.rad2deg(ishwr.zenith_true),np.rad2deg(ishwr.zenith_reco),np.rad2deg(ishwr.zenith_lapu))
    #   print("azimuth",np.rad2deg(ishwr.azimuth_true),np.rad2deg(ishwr.azimuth_reco),np.rad2deg(ishwr.azimuth_lapu))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_reco,ishwr.azimuth_reco)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_reco,ishwr.azimuth_reco,ishwr.zenith_lapu,ishwr.azimuth_lapu)))

    if sin_zenith == True:
      zenithAll = [np.sin(ishwr.zenith_true)**2 for ishwr in showerObjList]
      zenithAll_reco = [np.sin(ishwr.zenith_reco)**2 for ishwr in showerObjList]
      zenithAll_lapu = [np.sin(ishwr.zenith_lapu)**2 for ishwr in showerObjList]
      bins = np.linspace(0,1,21)
      xlabel = r"sin$^2$$\theta$"
    else:
      zenithAll = [np.rad2deg(ishwr.zenith_true) for ishwr in showerObjList]
      zenithAll_reco = [np.rad2deg(ishwr.zenith_reco) for ishwr in showerObjList]
      zenithAll_lapu = [np.rad2deg(ishwr.zenith_lapu) for ishwr in showerObjList]
      bins = np.linspace(0,90,91)
      xlabel = r"$\theta$$^{\circ}$"
    filtMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)>=zenCut]
    filtMisMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)<zenCut]
    impurity = sum(filtMisMatch)/(sum(filtMisMatch)+sum(filtMatch))
    weightAll = [ishwr.H4aWeight for ishwr in showerObjList]
    ax.hist(zenithAll,bins=bins,weights=weightAll,histtype="step",label="MC true")
    ax.hist(zenithAll_reco,bins=bins,weights=weightAll,histtype="step",label="plane")
    ax.hist(zenithAll_lapu,bins=bins,weights=weightAll,histtype="step",label="Laputop")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=18)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.text(0.06,0.95,s=r"impurity: {0:.1f}%".format(impurity*100),size=13,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
    ax.set_ylabel("count", fontsize=18)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_yscale("log")
    ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightTest.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightTest.png",transparent=False,bbox_inches='tight')
    #####################################################################
    ####################################################################
    # bins = np.linspace(0,90,91)
    # # bins = np.linspace(0,5,21)
    # openAngle_reco = [np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_reco,ishwr.azimuth_reco)) for ishwr in showerObjList]
    # openAngle_lapu = [np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_lapu,ishwr.azimuth_lapu)) for ishwr in showerObjList]
    # openAngle_lapuReco = [np.rad2deg(openingAngle(ishwr.zenith_reco,ishwr.azimuth_reco,ishwr.zenith_lapu,ishwr.azimuth_lapu)) for ishwr in showerObjList]
    # weightAll = [ishwr.H4aWeight for ishwr in showerObjList]
    # ax.hist(openAngle_reco,bins=bins,weights=weightAll,histtype="step",label="plane vs truth")
    # ax.hist(openAngle_lapu,bins=bins,weights=weightAll,histtype="step",label="Laputop vs truth")
    # ax.hist(openAngle_lapuReco,bins=bins,weights=weightAll,histtype="step",label="plane vs Laputop")
    # ax.tick_params(axis='both',which='both', direction='in', labelsize=11)
    # # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    # ax.set_ylabel("count", fontsize=18)
    # # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    # ax.set_xlabel(r"$\psi$ [$^{\circ}$]", fontsize=10)
    # ax.set_yscale("log")
    # ax.legend(fontsize="11")
    # plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightOpenTest.pdf",transparent=False,bbox_inches='tight')
    # plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightOpenTest.png",transparent=False,bbox_inches='tight')





  def plotHistWeightedFilterReco(self,triggerType,filterType):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(0,1,21)
    # showerObjList = [ishwr for ishwr in self.showerObjList if ishwr.nHLC < 1 and (ishwr.xcore**2+ishwr.ycore**2)<= 410**2]
    showerObjList = [ishwr for ishwr in self.showerObjList]
    zenithAll = [np.sin(ishwr.zenith_true)**2 for ishwr in showerObjList]
    weightAll = [ishwr.H4aWeight for ishwr in showerObjList]
    zenithAllWReco = [np.sin(ishwr.zenith_true)**2 for ishwr in showerObjList if not np.isnan(ishwr.zenith_reco) and not np.isnan(ishwr.zenith_lapu)]
    weightAllWReco = [ishwr.H4aWeight for ishwr in showerObjList if not np.isnan(ishwr.zenith_reco) and not np.isnan(ishwr.zenith_lapu)]
    zenithAll_reco = [np.sin(ishwr.zenith_reco)**2 for ishwr in showerObjList]
    weightAll_reco = [ishwr.H4aWeight for ishwr in showerObjList]
    zenithAll_reco45 = [np.sin(ishwr.zenith_reco)**2 for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45]
    weightAll_reco45 = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45]
    zenithAll_reco60 = [np.sin(ishwr.zenith_reco)**2 for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= 60]
    weightAll_reco60 = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= 60]
    zenithAll_lapu = [np.sin(ishwr.zenith_lapu)**2 for ishwr in showerObjList]
    weightAll_lapu = [ishwr.H4aWeight for ishwr in showerObjList]
    zenithAll_lapu45 = [np.sin(ishwr.zenith_lapu)**2 for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45]
    weightAll_lapu45 = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45]
    zenithAll_lapu60 = [np.sin(ishwr.zenith_lapu)**2 for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= 60]
    weightAll_lapu60 = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= 60]
    # ax.hist(zenithAll,bins=bins,weights=weightAll,histtype="step",label="MC true")
    ax.hist(zenithAllWReco,bins=bins,weights=weightAllWReco,histtype="step",label="MC true having both reco",alpha=0.6)
    ax.hist(zenithAll_reco,bins=bins,weights=weightAll_reco,histtype="step",label="plane",alpha=0.6)
    ax.hist(zenithAll_reco45,bins=bins,weights=weightAll_reco45,histtype="step",label="plane 45",alpha=0.6)
    ax.hist(zenithAll_reco60,bins=bins,weights=weightAll_reco60,histtype="step",label="plane 60",alpha=0.6)
    ax.hist(zenithAll_lapu,bins=bins,weights=weightAll_lapu,histtype="step",label="Laputop",alpha=0.6)
    ax.hist(zenithAll_lapu45,bins=bins,weights=weightAll_lapu45,histtype="step",label="Laputop 45",alpha=0.6)
    ax.hist(zenithAll_lapu60,bins=bins,weights=weightAll_lapu60,histtype="step",label="Laputop 60",alpha=0.6)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=11)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("count", fontsize=18)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel(r"sin$^2$$\theta$", fontsize=18)
    ax.set_yscale("log")
    ax.legend(fontsize="8")
    plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightFilt.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightFilt.png",transparent=False,bbox_inches='tight')

  def plotHistOpenAngle(self,triggerType,filterType,showerObjList):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(0,5,21)
    bins = np.linspace(-90,90,181)
    bins = np.linspace(0,20,41)
    zenCut = 45 #in degrees
    containment = True    
    sin_zenith = False
    if containment == True:
      showerObjList = [ishwr for ishwr in showerObjList if (ishwr.xcore**2+ishwr.ycore**2)<= 410**2]
    # showerObjList = [ishwr for ishwr in self.showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45 and np.rad2deg(ishwr.zenith_lapu) < 45]
    # showerObjList = [ishwr for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= zenCut and np.rad2deg(ishwr.zenith_true) < zenCut]
    showerObjList = [ishwr for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= zenCut]
    openAngle_reco = [np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_reco,ishwr.azimuth_reco)) for ishwr in showerObjList]
    openAngle_lapu = [np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_lapu,ishwr.azimuth_lapu)) for ishwr in showerObjList]
    openAngle_lapuReco = [np.rad2deg(openingAngle(ishwr.zenith_reco,ishwr.azimuth_reco,ishwr.zenith_lapu,ishwr.azimuth_lapu)) for ishwr in showerObjList]
    weightAll = [ishwr.H4aWeight for ishwr in showerObjList]
    ax.hist(openAngle_reco,bins=bins,weights=weightAll,histtype="step",label="plane vs truth")
    # ax.hist(openAngle_reco,weights=weightAll,histtype="step",label="plane vs truth")
    # ax.hist(openAngle_lapu,bins=bins,weights=weightAll,histtype="step",label="Laputop vs truth")
    # ax.hist(openAngle_lapuReco,bins=bins,weights=weightAll,histtype="step",label="plane vs Laputop")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=11)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("count", fontsize=18)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel(r"$\psi$ [$^{\circ}$]", fontsize=10)
    ax.set_yscale("log")
    ax.legend(fontsize="11")
    plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightOpen.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightOpen.png",transparent=False,bbox_inches='tight')

  def plotHistEnergy(self,triggerType,filterType):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(14,17,31)
    # energyAll_lowE = [ishwr.energy for ishwr in self.showerObjList if ishwr.energy*10**9<10**16]
    # energyAll_highE = [ishwr.energy for ishwr in self.showerObjList if ishwr.energy*10**9>=10**16]
    # energyTrigger_lowE = [ishwr.energy for ishwr in self.showerObjList if abs(getattr(ishwr,triggerType).value-1)<0.01 and ishwr.energy*10**9<10**16]
    # energyTrigger_highE = [ishwr.energy for ishwr in self.showerObjList if abs(getattr(ishwr,triggerType).value-1)<0.01 and ishwr.energy*10**9>=10**16]
    energyFilter = [np.log10(ishwr.energy*10**9) for ishwr in self.showerObjList if abs(getattr(ishwr,triggerType).value-1)<0.01 and abs(getattr(ishwr,filterType).value-1)<0.01]
    ax.hist(energyFilter,bins=bins,histtype="step",label="IceTopIncl_24 filt,E>=10PeV")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=11)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("count", fontsize=18)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel(r"lgE$", fontsize=18)
    # ax.set_yscale("log")
    ax.legend(fontsize="11")
    plt.savefig(plotFolder+"/hist"+str(filterType)+"Energy.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/hist"+str(filterType)+"Energy.png",transparent=False,bbox_inches='tight')


  def plotHistWeightedMisReco(self,triggerType,filterType,containment,showerObjList):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(0,1,21)
    zenCut = 45 #in degrees    
    sin_zenith = False
    if containment == True:
      showerObjList = [ishwr for ishwr in showerObjList if (ishwr.xcore**2+ishwr.ycore**2)<= 410**2]
    # showerObjList = [ishwr for ishwr in self.showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45 and np.rad2deg(ishwr.zenith_lapu) < 45]
    showerObjList = [ishwr for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= zenCut and np.rad2deg(ishwr.zenith_true) < zenCut]
    # for ishwr in showerObjList:
    #   print("zenith",np.rad2deg(ishwr.zenith_true),np.rad2deg(ishwr.zenith_reco),np.rad2deg(ishwr.zenith_lapu))
    #   print("azimuth",np.rad2deg(ishwr.azimuth_true),np.rad2deg(ishwr.azimuth_reco),np.rad2deg(ishwr.azimuth_lapu))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_reco,ishwr.azimuth_reco)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_reco,ishwr.azimuth_reco,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    if sin_zenith == True:
      zenithAll = [np.sin(ishwr.zenith_true)**2 for ishwr in showerObjList]
      zenithAll_reco = [np.sin(ishwr.zenith_reco)**2 for ishwr in showerObjList]
      zenithAll_lapu = [np.sin(ishwr.zenith_lapu)**2 for ishwr in showerObjList]
      bins = np.linspace(0,1,21)
      xlabel = r"sin$^2$$\theta$"
    else:
      zenithAll = [np.rad2deg(ishwr.zenith_true) for ishwr in showerObjList]
      zenithAll_reco = [np.rad2deg(ishwr.zenith_reco) for ishwr in showerObjList]
      zenithAll_lapu = [np.rad2deg(ishwr.zenith_lapu) for ishwr in showerObjList]
      bins = np.linspace(0,90,91)
      xlabel = r"$\theta$$^{\circ}$"
    filtMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)>=zenCut]
    filtMisMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)<zenCut]
    impurity = sum(filtMisMatch)/(sum(filtMisMatch)+sum(filtMatch))
    weightAll = [ishwr.H4aWeight for ishwr in showerObjList]
    ax.hist(zenithAll,bins=bins,weights=weightAll,histtype="step",label="MC true")
    ax.hist(zenithAll_reco,bins=bins,weights=weightAll,histtype="step",label="plane")
    # ax.hist(zenithAll_lapu,bins=bins,weights=weightAll,histtype="step",label="Laputop")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=18)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("count", fontsize=18)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.text(0.06,0.95,s=r"impurity: {0:.1f}%".format(impurity*100),size=13,horizontalalignment='left',verticalalignment='center', transform=ax.transAxes)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_yscale("log")
    ax.legend(fontsize="11")
    plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightMisReco.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/hist"+str(filterType)+"WeightMisReco.png",transparent=False,bbox_inches='tight')

  def plotHistWeightedMisRecoCores(self,triggerType,filterType,containment,showerObjList):    
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(0,1,21)
    zenCut = 45 #in degrees
    sin_zenith = False
    if containment == True:
      showerObjList = [ishwr for ishwr in showerObjList if (ishwr.xcore**2+ishwr.ycore**2)<= 410**2]
    # showerObjList = [ishwr for ishwr in self.showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45 and np.rad2deg(ishwr.zenith_lapu) < 45]
    showerObjList = [ishwr for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= zenCut and np.rad2deg(ishwr.zenith_true) < zenCut]
    # for ishwr in showerObjList:
    #   print("zenith",np.rad2deg(ishwr.zenith_true),np.rad2deg(ishwr.zenith_reco),np.rad2deg(ishwr.zenith_lapu))
    #   print("azimuth",np.rad2deg(ishwr.azimuth_true),np.rad2deg(ishwr.azimuth_reco),np.rad2deg(ishwr.azimuth_lapu))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_reco,ishwr.azimuth_reco)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_reco,ishwr.azimuth_reco,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    if sin_zenith == True:
      zenithAll = [np.sin(ishwr.zenith_true)**2 for ishwr in showerObjList]
      zenithAll_reco = [np.sin(ishwr.zenith_reco)**2 for ishwr in showerObjList]
      zenithAll_lapu = [np.sin(ishwr.zenith_lapu)**2 for ishwr in showerObjList]
      bins = np.linspace(0,1,21)
      xlabel = r"sin$^2$$\theta$"
    else:
      zenithAll = [np.rad2deg(ishwr.zenith_true) for ishwr in showerObjList]
      zenithAll_reco = [np.rad2deg(ishwr.zenith_reco) for ishwr in showerObjList]
      zenithAll_lapu = [np.rad2deg(ishwr.zenith_lapu) for ishwr in showerObjList]
      bins = np.linspace(0,90,91)
      xlabel = r"$\theta$$^{\circ}$"
    filtMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)>=zenCut]
    filtMisMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)<zenCut]
    impurity = sum(filtMisMatch)/(sum(filtMisMatch)+sum(filtMatch))
    core_X = [ishwr.xcore for ishwr in showerObjList]
    core_Y = [ishwr.ycore for ishwr in showerObjList]
    core_XCOG = [ishwr.xcoreCOG for ishwr in showerObjList]
    core_YCOG = [ishwr.ycoreCOG for ishwr in showerObjList]
    ax.plot(core_X,core_Y,"o",label="MC true")
    ax.plot(core_XCOG,core_YCOG,"d",label="COG")
    ax.set_aspect("equal")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=18)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("y [m]", fontsize=18)
    # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_xlabel(r"x [m]", fontsize=18)
    # ax.set_yscale("log")
    ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/"+str(filterType)+"MisRecoCore.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/"+str(filterType)+"MisRecoCore.png",transparent=False,bbox_inches='tight')

  def plotHistWeightedMisRecoQtot(self,triggerType,filterType,containment,showerObjList):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(0,1,21)
    zenCut = 45 #in degrees    
    sin_zenith = False
    if containment == True:
      showerObjList = [ishwr for ishwr in showerObjList if (ishwr.xcore**2+ishwr.ycore**2)<= 410**2]
    # showerObjList = [ishwr for ishwr in self.showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45 and np.rad2deg(ishwr.zenith_lapu) < 45]
    showerObjList = [ishwr for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= zenCut and np.rad2deg(ishwr.zenith_true) < zenCut]
    # for ishwr in showerObjList:
    #   print("zenith",np.rad2deg(ishwr.zenith_true),np.rad2deg(ishwr.zenith_reco),np.rad2deg(ishwr.zenith_lapu))
    #   print("azimuth",np.rad2deg(ishwr.azimuth_true),np.rad2deg(ishwr.azimuth_reco),np.rad2deg(ishwr.azimuth_lapu))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_reco,ishwr.azimuth_reco)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_reco,ishwr.azimuth_reco,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    if sin_zenith == True:
      zenithAll = [np.sin(ishwr.zenith_true)**2 for ishwr in showerObjList]
      zenithAll_reco = [np.sin(ishwr.zenith_reco)**2 for ishwr in showerObjList]
      zenithAll_lapu = [np.sin(ishwr.zenith_lapu)**2 for ishwr in showerObjList]
      bins = np.linspace(0,1,21)
      xlabel = r"sin$^2$$\theta$"
    else:
      zenithAll = [np.rad2deg(ishwr.zenith_true) for ishwr in showerObjList]
      zenithAll_reco = [np.rad2deg(ishwr.zenith_reco) for ishwr in showerObjList]
      zenithAll_lapu = [np.rad2deg(ishwr.zenith_lapu) for ishwr in showerObjList]
      bins = np.linspace(0,90,91)
      xlabel = r"$\theta$$^{\circ}$"
    # filtMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)>=zenCut]
    # filtMisMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)<zenCut]
    # impurity = sum(filtMisMatch)/(sum(filtMisMatch)+sum(filtMatch))
    Qtots_HLC = [ishwr.Qtot_HLC for ishwr in showerObjList]
    Qtots_SLC = [ishwr.Qtot_SLC for ishwr in showerObjList]
    Qtots_total = [ishwr.Qtot_SLC + ishwr.Qtot_HLC for ishwr in showerObjList]
    zenithAll = [np.rad2deg(ishwr.zenith_true) for ishwr in showerObjList]
    ax.plot(zenithAll,Qtots_HLC,"o",label="HLC")
    ax.plot(zenithAll,Qtots_SLC,"o",label="SLC")
    ax.plot(zenithAll,Qtots_total,"o",label="Sum")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=11)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    ax.set_xlabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("Q$_{tot} [VEM]$", fontsize=22)
    ax.set_yscale("log")
    ax.legend(fontsize="11")
    plt.savefig(plotFolder+"/"+str(filterType)+"MisRecoQtot.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/"+str(filterType)+"MisRecoQtot.png",transparent=False,bbox_inches='tight')


  def plotHistWeightedMisReconHit(self,triggerType,filterType,containment,showerObjList):    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(0,1,21)
    zenCut = 45 #in degrees    
    sin_zenith = False
    if containment == True:
      showerObjList = [ishwr for ishwr in showerObjList if (ishwr.xcore**2+ishwr.ycore**2)<= 410**2]
    # showerObjList = [ishwr for ishwr in self.showerObjList if np.rad2deg(ishwr.zenith_reco) >= 45 and np.rad2deg(ishwr.zenith_lapu) < 45]
    showerObjList = [ishwr for ishwr in showerObjList if np.rad2deg(ishwr.zenith_reco) >= zenCut and np.rad2deg(ishwr.zenith_true) < zenCut]
    # for ishwr in showerObjList:
    #   print("zenith",np.rad2deg(ishwr.zenith_true),np.rad2deg(ishwr.zenith_reco),np.rad2deg(ishwr.zenith_lapu))
    #   print("azimuth",np.rad2deg(ishwr.azimuth_true),np.rad2deg(ishwr.azimuth_reco),np.rad2deg(ishwr.azimuth_lapu))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_reco,ishwr.azimuth_reco)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_true,ishwr.azimuth_true,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    #   print(np.rad2deg(openingAngle(ishwr.zenith_reco,ishwr.azimuth_reco,ishwr.zenith_lapu,ishwr.azimuth_lapu)))
    if sin_zenith == True:
      zenithAll = [np.sin(ishwr.zenith_true)**2 for ishwr in showerObjList]
      zenithAll_reco = [np.sin(ishwr.zenith_reco)**2 for ishwr in showerObjList]
      zenithAll_lapu = [np.sin(ishwr.zenith_lapu)**2 for ishwr in showerObjList]
      bins = np.linspace(0,1,21)
      xlabel = r"sin$^2$$\theta$"
    else:
      zenithAll = [np.rad2deg(ishwr.zenith_true) for ishwr in showerObjList]
      zenithAll_reco = [np.rad2deg(ishwr.zenith_reco) for ishwr in showerObjList]
      zenithAll_lapu = [np.rad2deg(ishwr.zenith_lapu) for ishwr in showerObjList]
      bins = np.linspace(0,90,91)
      xlabel = r"$\theta$$^{\circ}$"
    # filtMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)>=zenCut]
    # filtMisMatch = [ishwr.H4aWeight for ishwr in showerObjList if np.rad2deg(ishwr.zenith_true)<zenCut]
    # impurity = sum(filtMisMatch)/(sum(filtMisMatch)+sum(filtMatch))
    Qtots_HLC = [ishwr.nHLC for ishwr in showerObjList]
    Qtots_SLC = [ishwr.nSLC for ishwr in showerObjList]
    Qtots_total = [ishwr.nHLC+ishwr.nSLC for ishwr in showerObjList]
    zenithAll = [np.rad2deg(ishwr.zenith_true) for ishwr in showerObjList]
    ax.plot(zenithAll,Qtots_HLC,"o",label="HLC")
    ax.plot(zenithAll,Qtots_SLC,"o",label="SLC")
    ax.plot(zenithAll,Qtots_total,"o",label="total")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=11)
    # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
    ax.set_xlabel(r"$\theta$$^{\circ}$", fontsize=22)
    ax.set_ylabel("nHit", fontsize=22)
    # ax.set_yscale("log")
    ax.legend(fontsize="11")
    plt.savefig(plotFolder+"/"+str(filterType)+"MisReconHit.pdf",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/"+str(filterType)+"MisReconHit.png",transparent=False,bbox_inches='tight')


  def Finish(self):
    # trigger = "IceTop7HGOnly"
    trigger = "IceTop7HGOnly"
    # trigger = "HLC6"
    if trigger == "IceTop7HGOnly":
      showerObjList = [ishwr for ishwr in self.showerObjList if abs(ishwr.HG7_3.value -1.0) < 0.1 and abs(ishwr.HLC6.value -0.0) < 0.1]
    elif trigger == "IceTop7HG":
      showerObjList = [ishwr for ishwr in self.showerObjList if abs(ishwr.HG7_3.value -1.0) < 0.1]
    elif trigger == "HLC6":
      showerObjList = [ishwr for ishwr in self.showerObjList if abs(ishwr.HLC6.value -1.0) < 0.1]
    # self.plotTrigEfficiency("HG7_3",containment=True)
    # self.plotFilterEfficiency("IceTopIncl_24",containment=True)
    # self.plotTrigEfficiencyZenith(triggerType="HG7_3",containment=True,wilson=True,weighting=True)
    # self.plotFilterEfficiencyZenith(filterType="IceTopIncl_24",containment=True,wilson=True,weighting=True)
    # self.plotEnergyFlux(triggerType="HG7_3",yscale="log",suffix="fluxLog",energyScale=0.0,containment=False)
    # self.plotZenithFlux("HG7_3","log","fluxLog",containment=False)
    # self.plotEnergyFluxFilter(triggerType="HG7_3",filterType="IceTopIncl_24",yscale="log",suffix="fluxLog",energyScale=0.0,containment=False)
    # self.plotZenithFluxFilter("HG7_3","IceTopIncl_24","log","fluxLog",containment=False)
    # self.plotHist("HG7_3","IceTopIncl_24")
    # self.plotHistWeighted("HG7_3","IceTopIncl_24")
    ###################################################################
    ###################################################################
    # self.plotHistWeightedTest("HG7_3","IceTopIncl_24",containment=True,showerObjList=showerObjList)
    # self.plotHistWeightedMisReco("HG7_3","IceTopIncl_24",containment=True,showerObjList=showerObjList)
    # self.plotHistWeightedMisRecoCores("HG7_3","IceTopIncl_24",containment=True,showerObjList=showerObjList)
    # self.plotHistWeightedMisRecoQtot("HG7_3","IceTopIncl_24",containment=True,showerObjList=showerObjList)
    # self.plotHistWeightedMisReconHit("HG7_3","IceTopIncl_24",containment=True,showerObjList=showerObjList)
    # self.plotHistWeightedFilterReco("HG7_3","IceTopIncl_24")
    self.plotHistOpenAngle("HG7_3","IceTopIncl_24",showerObjList)
    #####################################################
    #####################################################

def selectMisReco(frame,zenCut1,zenCut2):
  if abs(frame["HG7_3"].value -1.0) < 0.1 and abs(frame["HLC6"].value -0.0) < 0.1 and np.rad2deg(frame["ShowerPlaneIncl"].dir.zenith) > zenCut1 and np.rad2deg(frame["MCPrimary"].dir.zenith) < zenCut2 and (frame["MCPrimary"].pos.x**2+frame["MCPrimary"].pos.x**2)<= 410**2:
    return True
  else:
    return False



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

tray.AddModule(selectMisReco,"misReco",
              zenCut1 = 60,
              zenCut2 = 20,
              streams=[icetray.I3Frame.Physics],
              )

tray.AddModule("I3Writer","i3writer",
            filename=str(outputDir)+fileName,
            # Streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            DropOrphanStreams=[icetray.I3Frame.DAQ,icetray.I3Frame.Simulation],
            )

tray.Execute()
tray.Finish()