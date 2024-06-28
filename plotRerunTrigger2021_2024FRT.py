#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob

import datetime

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

from customColors import qualitative_colors

colorsCustom = qualitative_colors(5)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})


# #for 2023 January runs ( January 13 and January 21)
runs2024 = [int(irun) for irun in np.linspace(138858,138882,25)]
runs2023 = [int(irun) for irun in np.linspace(137540,137570,31) if int(irun) not in [137562,137555,137556,137557]]
runs2022 = [int(irun) for irun in np.linspace(136163,136190,28) if int(irun) not in [136167]]
runs2021 = [int(irun) for irun in np.linspace(134888,134915,28) if int(irun) not in [134892]]

runDict = {}
for irunName, irunList in zip(["2021","2022","2023","2024"],[runs2021,runs2022,runs2023,runs2024]):
  runDict[irunName] = irunList



runs = np.linspace(138858,138882,25)
# runs = [138867]
runs = [int(irun) for irun in runs]
####################################for 2024 frt rates##############################################
def getRate(year):
  rateFiles = sorted(glob.glob("/home/enpaudel/icecube/triggerStudy/FRTRate2024/triggerRateFRT2024??????.txt"))
  run = []
  HLCrate = []
  HG7rate = []
  for ifile in rateFiles:  
    with open(ifile,"r") as f:
      for line in f:
        splits = line.split(" ")
        run.append(int(splits[0]))
        HG7rate.append(int(splits[1])*100.0/int(splits[3]))
        HLCrate.append(int(splits[2])*100.0/int(splits[3]))
  return HG7rate
################for 2023 December ########################################





fig = plt.figure(figsize=(8,5))
gs = gridspec.GridSpec(nrows=1,ncols=1)
ax = fig.add_subplot(gs[0])
# ax.errorbar(1,np.mean(HG7rateList),yerr=np.std(HG7rateList)/np.sqrt(len(HG7rateList)),fmt="s",ls="-",lw = 2.5,ms=10,c=colorsCustom[2],label=r"{}".format("data"),alpha=0.8)
ax.errorbar(3,np.mean(HG7rate),yerr=np.std(HG7rate)/np.sqrt(len(HG7rate)),fmt="o",ls="-",lw = 2.5,ms=10,c=colorsCustom[0],label=r"{}".format(""),alpha=0.8)
ax.tick_params(axis='both',which='both', direction='in', labelsize=18)
ax.set_ylabel("rate [Hz]", fontsize=20)
ax.set_xlabel("time", fontsize=20)
ax.set_xticks([1,3],["2023 December","2024 January"])
ax.legend(fontsize=16,ncol=2)
ax.set_ylim(0,35)
ax.set_xlim(0,4)
ax.grid(True,alpha=0.5)
plt.savefig(plotFolder+"/triggerRateFRT2021_2024.pdf",transparent=False,bbox_inches='tight')
plt.savefig(plotFolder+"/triggerRateFRT2021_2024.png",transparent=False,bbox_inches='tight')