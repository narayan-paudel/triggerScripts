#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube.dataclasses import I3Constants

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)


import scipy.interpolate



import glob

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--vertEvents", dest="vertEvents",
                    default=False, action="store_true", required=False,
                    help="plot vertical events")
args = parser.parse_args()

# usage python hitsInclinedShower.py --vertEvents


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

inFile = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz"
# fileList = sorted(glob.glob(fileDir+"*.i3.*"))[:2]
# fileDir = "/data/user/enpaudel/triggerStudy/simFiles/dataSetClean1_6/"
fileList = sorted(glob.glob(fileDir+"*.i3.*"))[:4]

inclinationCut = 60 #degree
energyCut = 10**16 #eV

exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

# outputDir = "/home/enpaudel/dataExp/dataSetClean_InclinedHE/"

outputDir = "/home/enpaudel/dataExp/dataSetClean_InclinedHE_7HG_reco"
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

def openingAngle(theta1,phi1,theta2,phi2):
  return np.arccos(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2))

def test7HG(frame):
  return frame["HG7_3"]>0

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



# def getSPTime(particle,omkey,pulseTime):
#   '''
#   returns a shower plane time. Both script works exactly same
#   '''
#   domgeo = geoFrame["I3Geometry"].omgeo[omkey]
#   x_om = domgeo.position.x
#   y_om = domgeo.position.y
#   z_om = domgeo.position.z
#   t_offset = (particle.dir.x * (x_om - particle.pos.x) + \
#                                   particle.dir.y * (y_om - particle.pos.y) + \
#                                   particle.dir.z * (z_om - particle.pos.z))/dataclasses.I3Constants.c
#   print("t_offset",t_offset+pulseTime)
#   return pulseTime + t_offset



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

class hitObj(object):
  """docstring for hitObj"""
  def __init__(self, time,x,y,z,charge):
    super(hitObj, self).__init__()
    self.time = time
    self.x = x
    self.y = y
    self.z = z
    self.charge = charge
    # self.dplane = self.distanceToPlane()
  def shiftCOG(self,xCOG,yCOG,zCOG,tCOG):  
    self.x -= xCOG
    self.y -= yCOG
    self.z -= zCOG
    self.time -= tCOG

  def distanceToPlane(self,nx,ny,nz):
    '''
    when the plane pass through the origin
    point(x,y,z) to plane (nx,ny,nz)
    '''
    # return (x*nx+y*ny+z*nz)/np.sqrt(nx**2+ny**2+ny**2)
    self.distPlane = (self.x*nx+self.y*ny+self.z*nz)/np.sqrt(nx**2+ny**2+ny**2)


class PlaneFit(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)

  def Configure(self):
    self.ZenReco = []
    self.count = 0
    self.n_cog = 7
    self.power = 1

  def Geometry(self,frame):
    self.domgeo = frame["I3Geometry"].omgeo

  def getHitCOG(self,hitList,n_cog):
    if len(hitList)>n_cog:
      coghitList = hitList[:self.n_cog]
    else:
      coghitList = hitList
    weightTotal = np.sum([ihit.charge**self.power for ihit in coghitList])
    meanX = np.sum([ihit.x * ihit.charge**self.power for ihit in coghitList])/weightTotal
    meanY = np.sum([ihit.y * ihit.charge**self.power for ihit in coghitList])/weightTotal
    meanZ = np.sum([ihit.z * ihit.charge**self.power for ihit in coghitList])/weightTotal
    meant = np.sum([ihit.time * ihit.charge**self.power for ihit in coghitList])/weightTotal
    meanCharge = np.sum([ihit.charge * ihit.charge**self.power for ihit in coghitList])/weightTotal
    meanHit = hitObj(meant,meanX,meanY,meanZ,meanCharge)
    return meanHit
  def getUnitPlane(self,hitList):
    '''
    returns nx, ny, nz, unit vector for plane
    '''
    x = [ihit.x for ihit in hitList]
    y = [ihit.y for ihit in hitList]
    z = [ihit.z for ihit in hitList]
    time = [ihit.time for ihit in hitList]
    Sxx = np.sum([ix*ix for ix in x])
    Syy = np.sum([iy*iy for iy in y])
    Szz = np.sum([iz*iz for iz in z])
    Sxy = np.sum([ix*iy for ix,iy in zip(x,y) ])
    Szx = np.sum([ix*iz for ix,iz in zip(x,z) ])
    Syz = np.sum([iy*iz for iy,iz in zip(y,z) ])
    Stx = np.sum([ix*I3Constants.c*it/I3Units.second for ix,it in zip(x,time) ])
    Sty = np.sum([iy*I3Constants.c*it/I3Units.second for iy,it in zip(y,time) ])
    Stz = np.sum([iz*I3Constants.c*it/I3Units.second for iz,it in zip(z,time) ])
    S = np.array([[Sxx,Sxy,Szx],[Sxy,Syy,Syz],[Szx,Syz,Szz]])
    T = [[Stx],[Sty],[Stz]]
    detS = np.linalg.det(S)
    # print("matrix S",S)
    # print("det matrix S",detS)
    if abs(detS)>0:
      invS = np.linalg.inv(S)
      invST = np.matmul(invS,T)
      # print("detS",detS,invS,invST,invST.T)
      nx,ny,nz = invST.T[0]
      mag = np.sqrt(nx**2+ny**2+nz**2)
      return nx/mag,ny/mag,nz/mag      
    else:
      return np.nan,np.nan,np.nan

  def removeFurthestHit(self,hitList):
    distancesList = [abs(ihit.distPlane) for ihit in hitList]
    maxD = max(distancesList)
    updatedHitList = []
    for ihit in hitList:
      if abs(maxD - abs(ihit.distPlane)) > 0.05:
        updatedHitList.append(ihit)
    return updatedHitList


  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    # qtotHLC = frame["OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalCharge"]
    # qtotSLC = frame["OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalCharge"]
    # qtot = frame["IceTopVEMPulsesTotalCharge"]
    hitList = []
    psm = frame["IceTopVEMPulses"]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap or psm.__class__ == dataclasses.I3RecoPulseSeriesMapUnion:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, "IceTopVEMPulses")
    for om,pulses in psm:
      if om in hg_doms:
        x = self.domgeo[om].position.x
        y = self.domgeo[om].position.y
        z = self.domgeo[om].position.z
        if len(pulses)>1:
          print("multiple pulse per DOM",len(pulses))
        for pulse in pulses:
          time = pulse.time
          charge = pulse.charge
          break
        hitList.append(hitObj(time,x,y,z,charge))
    hitList = sorted(hitList, key=lambda x:x.charge,reverse=True)[:7]
    hitCOG = self.getHitCOG(hitList,self.n_cog)
    for ihit in hitList:
      ihit.shiftCOG(hitCOG.x,hitCOG.y,hitCOG.z,hitCOG.time)
    while len(hitList) > 3:
      print("iter",len(hitList))
      nx,ny,nz = self.getUnitPlane(hitList)
      if nx is not np.isnan and ny is not np.isnan and nz is not np.isnan:
        theta = np.rad2deg(np.arctan(ny/nx))
        phi = np.rad2deg(np.arccos(nx))
      else:
        theta = np.nan
        phi = np.nan
        break
      # phi = np.rad2deg(np.arccos(nx/(np.sqrt(nx**2+ny**2+nz**2))))
      for ihit in hitList:
        ihit.distanceToPlane(nx,ny,nz)
      hitList = self.removeFurthestHit(hitList)
      print("theta",np.rad2deg(frame["MCPrimary"].dir.zenith),np.rad2deg(frame["ShowerPlane"].dir.zenith),theta)
      print("phi",np.rad2deg(frame["MCPrimary"].dir.azimuth),np.rad2deg(frame["ShowerPlane"].dir.azimuth),phi)
    frame["MyRecoTheta"] = dataclasses.I3Double(np.deg2rad(theta))
    frame["MyRecoPhi"] = dataclasses.I3Double(np.deg2rad(phi))
    self.PushFrame(frame)

    self.count += 1
    # print("count",self.count)
    if self.count < 10:
      self.fig = plt.figure(figsize=(8,5))
      gs = gridspec.GridSpec(nrows=1,ncols=1)
      self.ax = self.fig.add_subplot(gs[0])
      x = [ihit.x for ihit in hitList]
      y = [ihit.y for ihit in hitList]
      z = [ihit.z for ihit in hitList]
      self.ax.plot(x,y,".",label=r"",alpha=1)
      self.ax.plot(hitCOG.x,hitCOG.y,"o")
      self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
      # self.ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
      # self.ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
      self.ax.set_ylabel("x", fontsize=22)
      # self.ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
      self.ax.set_xlabel("y", fontsize=22)
      # self.ax.set_yscale("log")
      self.ax.set_xlim(-600,600)
      self.ax.set_ylim(-600,600)
      self.ax.legend()
      plt.savefig(plotFolder+"/hits"+str(self.count)+".pdf",transparent=False,bbox_inches='tight')
      plt.savefig(plotFolder+"/hits"+str(self.count)+".png",transparent=False,bbox_inches='tight')



  def Finish(self):
    pass


class zenithCheck(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.openingAngleList = []
    self.openingAngleMyRecoList = []
    self.r_diffList = []


  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    xcore = frame["MCPrimary"].pos.x
    ycore = frame["MCPrimary"].pos.y
    xcore_reco = frame["ShowerCOG"].pos.x
    ycore_reco = frame["ShowerCOG"].pos.y
    r_diff = np.sqrt((xcore-xcore_reco)**2+(ycore-ycore_reco)**2)/I3Units.m
    # print("cores",(xcore,ycore),(xcore_reco,ycore_reco),(xcore-xcore_reco,ycore-ycore_reco),r_diff)
    zenith_true = frame["MCPrimary"].dir.zenith
    azimuth_true = frame["MCPrimary"].dir.azimuth
    zenith_reco = frame["ShowerPlane"].dir.zenith
    azimuth_reco = frame["ShowerPlane"].dir.azimuth
    openAngle = openingAngle(zenith_true,azimuth_true,zenith_reco,azimuth_reco)
    openAngleMyReco = openingAngle(zenith_true,azimuth_true,frame["MyRecoTheta"].value,frame["MyRecoPhi"].value)
    # print("zenith True",zenith_reco,zenith_true,np.arcsin(np.sqrt(self.zenithBin[0])),np.arcsin(np.sqrt(self.zenithBin[1])))
    # if np.arcsin(np.sqrt(self.zenithBin[0])) <= zenith_true < np.arcsin(np.sqrt(self.zenithBin[1])) and not np.isnan(openAngle):
    self.openingAngleList.append(openAngle*180.0/np.pi)
    self.openingAngleMyRecoList.append(openAngleMyReco*180.0/np.pi)
    self.r_diffList.append(r_diff)

  def Finish(self):
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    bins = np.linspace(-1,100,102)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    self.ax.hist(self.openingAngleList,bins=bins,histtype="step",label=r"SimplePlane",lw=2.5)
    self.ax.hist(self.openingAngleMyRecoList,bins=bins,histtype="step",label=r"MyReco",lw=2.5)
    self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\psi$ [$^{\circ}$]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_xlim(0,100)
    # self.ax.set_ylim(0,4300)
    self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/openAngleHG7OnlyReco2{}.png".format(plotSuffix),transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/openAngleHG7OnlyReco2{}.pdf".format(plotSuffix),transparent=False,bbox_inches='tight')
    plt.close()
    #####################################################
    #####################################################


tray = I3Tray()
tray.AddModule("I3Reader","reader",
             # filenameList = args.input,
             filenameList = [GCD]+fileList,
             # filename = inFile,
            )

tray.AddModule(test7HG,"7HG",
              streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
              )

def Unify(frame, Keys, Output):
  """
  Simple utility to merge RecoPulseSerieses into a single Union.
  """
  extants = [k for k in Keys if k in frame]
  union = dataclasses.I3RecoPulseSeriesMapUnion(frame, extants)
  frame[Output] = union
tray.Add(Unify,"UnionHLCSLC",
  # Keys=["OfflineIceTopHLCTankPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge"],
  Keys=["OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge"],
  # Output='IceTopTankPulses',
  Output='IceTopVEMPulses',
  streams = [icetray.I3Frame.DAQ],
  )
tray.AddModule(PlaneFit,"pFit",
               # Keys=["IceTopTankPulses"],
               # Keys=["IceTopVEMPulses"],
                )
tray.AddModule(zenithCheck,"zhits",
            # streams = [icetray.I3Frame.DAQ],
            )

tray.AddModule("I3Writer","i3writer",
            filename=str(outputDir)+fileName,
            streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            )

tray.Execute()
tray.Finish()