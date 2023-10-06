#!/usr/bin/env python3

from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('gcd', type=str, default="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz", help='Input GCD')
args = parser.parse_args()

original_elevation = 883.92
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

def snowheight(tank,snowheight):
  """additionally adds a circle of radius r m """
  fig = plt.figure(figsize=(40,8))
  gs = gridspec.GridSpec(ncols=1,nrows=1)
  ax = fig.add_subplot(gs[0])
  # ax.plot(tank,snowheight,"-",ms=7,mew=2,mfc='none',label="snow height",alpha=1)
  ax.plot(tank,snowheight,"-",label="snow height",alpha=1)
  ax.set_xlabel(r"tank", fontsize=20)
  ax.set_ylabel(r"snow height [m]", fontsize=20)
  ax.tick_params(axis='both',which='both',direction='in', labelsize=10)
  ax.tick_params(which='both', width=1.5)
  # ax.tick_params(which='major', length=7)
  # ax.tick_params(which='minor', length=4)
  # ax.set_yscale('log')
  ax.grid(True,alpha=0.5)
  ax.legend(fontsize=14)
  # ax.legend(fontsize=14,ncol=2)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  plt.savefig(plotFolder+"tankSnowHeight.pdf",transparent=False,bbox_inches='tight')
  plt.close()
snowHeight = []
for frame in dataio.I3File(args.gcd,'r'):
  if frame.Stop == icetray.I3Frame.Geometry:
    geom = frame['I3Geometry']
    stageo = geom.stationgeo
    x = []
    sh = []
    for e,st in stageo:
      for tank in st:
        for om in tank.omkey_list:
          if om[1] in [61,63]:
            x.append("{:d}\n{}".format(e,om[1]).replace("61","A").replace("63","B"))
            sh.append(original_elevation+tank.position.z+tank.snowheight)
            snowHeight.append(tank.snowheight)
print(max(snowHeight),min(snowHeight))

snowheight(x,sh)