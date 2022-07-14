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

from customColors import qualitative_colors

import numpy as np
from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
import os
import numpy
# gcd_file = "/home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz"
gcd_file = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz"

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"


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
          if omkey[1] == 61:
            x[0].append(tank.position.x)
            y[0].append(tank.position.y)
            z[0].append(tank.position.z)
          elif omkey[1] == 63:
            x[1].append(tank.position.x)
            y[1].append(tank.position.y)
            z[1].append(tank.position.z)

def getCircle(radius):
  theta = np.linspace(0,2*np.pi,100)
  return radius*np.cos(theta),radius*np.sin(theta)

def tankLayout(x,y,stations,r):
  """additionally adds a circle of radius r m """
  fig = plt.figure(figsize=(8,8))
  gs = gridspec.GridSpec(ncols=1,nrows=1)
  ax = fig.add_subplot(gs[0])
  ax.plot(x[0],y[0],"o",ms=7,mew=2,mfc='none',c=qualitative_colors(5)[1],label="tank A",alpha=1)
  ax.plot(x[1],y[1],"o",ms=7,mew=2,mfc='none',c=qualitative_colors(5)[2],label="tank B",alpha=1)
  xCirc,yCirc = getCircle(r)
  ax.plot(xCirc,yCirc,'-',c=qualitative_colors(5)[3],lw=3.0,label="r = {:.1f} m".format(r))
  ax.set_xlabel(r"x [m]", fontsize=20)
  ax.set_ylabel(r"y [m]", fontsize=20)
  ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
  ax.tick_params(which='both', width=1.5)
  ax.tick_params(which='major', length=7)
  ax.tick_params(which='minor', length=4)
  # ax.set_yscale('log')
  ax.grid(True,alpha=0.5)
  ax.set_aspect("equal")
  ax.set_ylim(-650,650)
  ax.set_xlim(-650,650)
  for ix,iy,ista in zip(x[0],y[0],stations):
    ax.text(ix-40,iy-40,s=r"{}".format(ista),size=12)
  ax.legend(fontsize=14)
  # ax.legend(fontsize=14,ncol=2)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  plt.savefig(plotFolder+"tankLayout{:d}.pdf".format(r),transparent=False,bbox_inches='tight')
  plt.close()

tankLayout(x,y,stations,640)

         