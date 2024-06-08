#!/usr/bin/env python3

import os
import tables
import matplotlib as mpl
mpl.use("Agg")
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
from matplotlib.ticker import AutoMinorLocator,MultipleLocator
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm

from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.markers import MarkerStyle
from matplotlib.transforms import Affine2D

from customColors import qualitative_colors

import numpy as np
from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
import os
import numpy
# gcd_file = "/home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz"
gcd_file = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

inputFile = "/home/enpaudel/icecube/triggerStudy/simFilesTest/inclinedFe.i3.gz"

eventID = 5421454
pulses = "IceTopPulses"

posKeyDict = {}




def semicircle_path(radius=1, direction='up'):
    # Define semicircle points
    theta = np.linspace(0, np.pi, 100)
    if direction == 'down':
        theta = np.linspace(np.pi, 2*np.pi, 100)
    elif direction == 'left':
        theta = np.linspace(0.5*np.pi, 1.5*np.pi, 100)
    elif direction == 'right':
        theta = np.linspace(-0.5*np.pi, 0.5*np.pi, 100)
        
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    # Close the path to make a filled marker
    if direction in ['up', 'down']:
        x = np.append(x, [0])
        y = np.append(y, [0])
    elif direction in ['left', 'right']:
        x = np.append(x, [0])
        y = np.append(y, [0])
        
    vertices = np.vstack((x, y)).T
    codes = [Path.MOVETO] + [Path.LINETO] * (len(vertices) - 1) + [Path.CLOSEPOLY]
    path = Path(vertices, codes)
    return path

# Create a custom marker
def semicircle_marker(direction='up'):
    semicircle = semicircle_path(direction=direction)
    return PathPatch(semicircle, transform=plt.gca().transData)


for f in dataio.I3File(gcd_file,'r'):
  if f.Stop == icetray.I3Frame.Geometry:
    geom = f['I3Geometry']
    # print("geometry",geom.keys())
    stageo = geom.stationgeo
    omgeo = geom.omgeo
    tankgeo = geom.tankgeo
    x = [[],[]]
    y = [[],[]]
    z = [[],[]]
    stations = []
    omkeys = []
    for stnkey,station in geom.stationgeo:
      stations.append(stnkey)
      # print("station key",stnkey,station)
      for tank in station:
        for omkey in tank.omkey_list:
          posKeyDict[omkey] = tank.position
          omkeys.append(omkey)


def getxy(om,posKeyDict):
  return posKeyDict[om]

def getPlotParams(inputFile,eventID,pulses):
  for f in dataio.I3File(inputFile,'r'):
    if f.Stop == icetray.I3Frame.DAQ:
      if f['I3EventHeader'].event_id == eventID:      
        eventN = f['I3EventHeader'].event_id
        runN = f['I3EventHeader'].run_id
        ps = f[pulses]
        psm = dataclasses.I3RecoPulseSeriesMap.from_frame(f,pulses)
        markerListA = []
        markerListB = []
        xAlist = []
        yAlist = []
        tAlist = []
        qAlist = []
        xBlist = []
        yBlist = []
        tBlist = []
        qBlist = []
        om_hits = []
        for om,pulses in psm:
          om_hits.append(om)
          if om.om in [61,62]:
            xA = posKeyDict[om].x
            yA = posKeyDict[om].y
            chargeA = sum([c.charge for c in pulses])
            tA = pulses[0].time*I3Units.ns/1000
            xAlist.append(xA)
            yAlist.append(yA)
            tAlist.append(tA)
            qAlist.append(np.log10(chargeA))
            markerListA.append(MarkerStyle("o", fillstyle="left"))
          elif om.om in [63,64]:
            xB = posKeyDict[om].x
            yB = posKeyDict[om].y
            chargeB = sum([c.charge for c in pulses])
            tB = pulses[0].time*I3Units.ns/1000
            if eventID == 741466:              
              if tB < 373.78:
                xBlist.append(xB)
                yBlist.append(yB)
                tBlist.append(tB)
                qBlist.append(np.log10(chargeB))
                markerListB.append(MarkerStyle("o", fillstyle="right"))
            else:
              xBlist.append(xB)
              yBlist.append(yB)
              tBlist.append(tB)
              qBlist.append(np.log10(chargeB))
              markerListB.append(MarkerStyle("o", fillstyle="right"))
        print("id",eventN,runN)
  print(om_hits)
  om_unhits = [om for om in omkeys if om not in om_hits]
  xA_unhit = [posKeyDict[om].x for om in om_unhits if om.om in [61]]
  xB_unhit = [posKeyDict[om].x for om in om_unhits if om.om in [63]]
  yA_unhit = [posKeyDict[om].y for om in om_unhits if om.om in [61]]
  yB_unhit = [posKeyDict[om].y for om in om_unhits if om.om in [63]]
  print(om_hits)
  print("maxt",max(tBlist))
  # tBlist.remove(max(tBlist))
  tlist = tAlist + tBlist
  print("maxt",max(tBlist),max(tAlist),max(tlist),min(tlist),max(tlist)-min(tlist))
  qlist = qAlist + qBlist
  tAlist_n = (tAlist - np.min(tlist)) / (np.max(tlist) - np.min(tlist))
  tBlist_n = (tBlist - np.min(tlist)) / (np.max(tlist) - np.min(tlist))
  qAlist_n = (qAlist - np.min(qlist)) / (np.max(qlist) - np.min(qlist))
  qBlist_n = (qBlist - np.min(qlist)) / (np.max(qlist) - np.min(qlist))
  return xAlist,yAlist,xBlist,yBlist,xA_unhit,yA_unhit,xB_unhit,yB_unhit,tlist,tAlist_n,tBlist_n,qAlist_n,qBlist_n,eventID

def eventView(xA,yA,xB,yB,xA_unhit,yA_unhit,xB_unhit,yB_unhit,tlist,tAlist,tBlist,qA,qB,eventID):
  """additionally adds a circle of radius r m """
  # colors = cmap(tlist)
  fig = plt.figure(figsize=(8,8))
  gs = gridspec.GridSpec(ncols=1,nrows=1)
  ax = fig.add_subplot(gs[0])
  # ax.plot(x,y,"o",ms=7,mew=2,mfc='none',c=cmap,s=q,label="tank A",alpha=1)
  # print("plotting scatter",len(x),len(y),len(colors),len(q))
  cmap = plt.get_cmap('rainbow_r')
  colorsA = cmap(tAlist_n)
  colorsB = cmap(tBlist_n)
  qA = np.asarray(qA)
  qB = np.asarray(qB)
  # marker = semicircle_marker('up')
  # ax.scatter(xA,yA,c=colors,s=1000*q, marker=MarkerStyle("o", fillstyle="left"),label="tank A",alpha=0.5)
  t = Affine2D().rotate_deg(45)
  sc=ax.scatter(xA,yA,c=colorsA,s=1500*qA, marker=MarkerStyle("o", fillstyle="left"),label="tank A",alpha=1)
  ax.scatter(xB,yB,c=colorsB,s=1500*qB, marker=MarkerStyle("o", fillstyle="right"),label="tank B",alpha=1)
  ax.scatter(xA_unhit,yA_unhit,marker=MarkerStyle("o", fillstyle="left"),c="gray",label="unhit",alpha=0.3)
  ax.scatter(xB_unhit,yB_unhit,marker=MarkerStyle("o", fillstyle="right"),c="gray",label="unhit",alpha=0.3)
  # ax.scatter(xB_unhit,yB_unhit,marker=MarkerStyle("o", fillstyle="right"),ms=5,mew=2,mfc='none',c="gray",label="unhit",alpha=0.3)
  # ax.plot(x[1],y[1],"o",ms=7,mew=2,mfc='none',c=qualitative_colors(5)[2],label="tank B",alpha=1)
  # ax.plot(xCirc,yCirc,'-',c=qualitative_colors(5)[3],lw=3.0,label="r = {:.1f} m".format(r))
  ax.set_xlabel(r"x [m]", fontsize=24)
  ax.set_ylabel(r"y [m]", fontsize=24)
  ax.tick_params(axis='both',which='both',direction='in', labelsize=24)
  ax.tick_params(which='both', width=1.5)
  ax.tick_params(which='major', length=7)
  ax.tick_params(which='minor', length=4)
  # ax.set_yscale('log')
  ax.grid(True,alpha=0.5)
  ax.set_aspect("equal")
  ax.set_ylim(-650,650)
  ax.set_xlim(-650,650)
  # for ix,iy,ista in zip(x[0],y[0],stations):
  #   ax.text(ix-40,iy-40,s=r"{}".format(ista),size=12)
  # ax.legend(fontsize=14)
  # ax.legend(fontsize=14,ncol=2)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  timeBar = [0.0, 0.5, 1.0]
  timeBarScaled = [ielt*(np.max(tlist) - np.min(tlist)) for ielt in timeBar]
  cbar = fig.colorbar(sc,ax=ax,ticks=timeBar,fraction=0.047,pad=0.015)
  sc.set_cmap('rainbow_r')
  print(timeBar,["{:.2f}".format(timeBarScaled[0]),"{:.2f}".format(timeBarScaled[1]),"{:.2f}".format(timeBarScaled[2])])
  # cbar.ax.set_xticklabels([ielt*(np.max(tlist) - np.min(tlist)) for ielt in timeBar])
  cbar.set_ticklabels(["{:.2f}".format(timeBarScaled[0]),"{:.2f}".format(timeBarScaled[1]),"{:.2f}".format(timeBarScaled[2])])
  cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
  cbar.set_label(r'time [$\mathrm{\mu}$s]',fontsize=24)
  # cbar.ax.set_ylabel(r"count", fontsize=24)
  plt.savefig(plotFolder+"eventView{:d}.pdf".format(eventID),transparent=False,bbox_inches='tight')
  plt.savefig(plotFolder+"eventView{:d}.png".format(eventID),transparent=False,bbox_inches='tight')
  plt.close()

inputFile = "/home/enpaudel/icecube/triggerStudy/simFilesTest/inclinedFe.i3.gz"
eventID = 5421454
xAlist,yAlist,xBlist,yBlist,xA_unhit,yA_unhit,xB_unhit,yB_unhit,tlist,tAlist_n,tBlist_n,qAlist_n,qBlist_n,eventID = getPlotParams(inputFile,eventID,pulses)
eventView(xAlist,yAlist,xBlist,yBlist,xA_unhit,yA_unhit,xB_unhit,yB_unhit,tlist,tAlist_n,tBlist_n,qAlist_n,qBlist_n,eventID)

inputFile = "/home/enpaudel/icecube/triggerStudy/simFilesTest/verticalFe.i3.gz"
eventID = 741466
xAlist,yAlist,xBlist,yBlist,xA_unhit,yA_unhit,xB_unhit,yB_unhit,tlist,tAlist_n,tBlist_n,qAlist_n,qBlist_n,eventID = getPlotParams(inputFile,eventID,pulses)
eventView(xAlist,yAlist,xBlist,yBlist,xA_unhit,yA_unhit,xB_unhit,yB_unhit,tlist,tAlist_n,tBlist_n,qAlist_n,qBlist_n,eventID)

     