#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob

import datetime

import pandas as pd

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
weekwise = True

class runObj(object):
  """docstring for runObj"""
  def __init__(self, runID, runDur):
    super(runObj, self).__init__()
    self.runID = runID
    self.runDur = runDur


if weekwise:
  ##########################################################################################################
  GRL = "/data/exp/IceCube/2023/filtered/level2/"
  GRL2023 = GRL + "IC86_2023_GoodRunInfo.txt"
  GRL2023_df = pd.read_csv(GRL2023,sep="\\s+",header=0,escapechar='#')
  GRL_IT_List = np.asarray(GRL2023_df.iloc[:,0][GRL2023_df.iloc[:,2]==1])
  liveTime = np.asarray(GRL2023_df.iloc[:,3][GRL2023_df.iloc[:,2]==1])
  runList = [int(irun) for irun in GRL_IT_List]
  runTimeDict = {}
  for irun,idur in zip(runList,liveTime):
    runTimeDict[irun] = idur
  removeList = [138629,138640,138657,138658,138659,138734,138764,138768,138676,138679,138696,
  138712,138744,138789,138823,138824,138838,138883]
  for ielt in removeList:
    runList.remove(ielt)
  runObjList = []
  ######################measured rate###############################
  for irun in runList:
    thisRun = runObj(irun,runTimeDict[irun])
    # with open("/home/enpaudel/icecube/triggerStudy/triggerRateListCombinedShort{}.txt".format(int(run)),"r") as f:
    # print("this run",irun)
    with open("/home/enpaudel/icecube/triggerStudy/rateFiles/triggerRate{}.txt".format(int(irun)),"r") as f:
      iHG7rate = []
      iHLCrate = []
      for line in f:
        splits = line.split(" ")
        if len(splits) > 2:
          iHG7rate.append(int(splits[2]))
          iHLCrate.append(int(splits[3]))
    # print("rates of {} {:.2f} Hz".format(irun,sum(iHG7rate)/runTimeDict[irun]))
    thisRun.HG7rate = sum(iHG7rate)/thisRun.runDur
    thisRun.HLCrate = sum(iHLCrate)/thisRun.runDur
    with open("/home/enpaudel/icecube/triggerStudy/FRTRate2023/triggerRateFRT2023_{}.txt".format(int(irun)),"r") as f:
      for line in f:
        splits = line.split(" ")
        # run.append(int(splits[0]))
    thisRun.HG7rateFRT = int(splits[1])*100.0/int(splits[3])
    thisRun.HLCrateFRT = int(splits[2])*100.0/int(splits[3])
    runObjList.append(thisRun)
  weekDict = {"12/03":138634,"12/10":138666,"12/17":138693,"12/24":138756,
  "12/31":138805,"01/07":138830,"01/14":138861,"01/21":138882}
  weekKeys = list(weekDict.keys())

  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])

  for i,iweek in enumerate(weekKeys[:-1]):
    startRun = weekDict[weekKeys[i]]
    endRun = weekDict[weekKeys[i+1]]
    runThisWeek = [irun for irun in runObjList if startRun <= irun.runID < endRun]
    HG7rateList = [irun.HG7rate for irun in runThisWeek]
    HG7rateFRTList = [irun.HG7rateFRT for irun in runThisWeek]
    ax.errorbar(iweek,np.mean(HG7rateList),yerr=np.std(HG7rateList)/np.sqrt(len(HG7rateList)),fmt="s",ls="-",lw = 2.5,ms=12,c=colorsCustom[0],label='_nolegend_' if i>0 else r"{}".format("data"),alpha=0.8)
    ax.errorbar(iweek,np.mean(HG7rateFRTList),yerr=np.std(HG7rateFRTList)/np.sqrt(len(HG7rateFRTList)),fmt="o",ls="-",lw = 2.5,ms=12,c=colorsCustom[2],label='_nolegend_' if i>0 else r"{}".format("FRT"),alpha=0.8)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=20, pad=8)
  ax.set_ylabel("rate [Hz]", fontsize=20)
  ax.set_xlabel("week", fontsize=20)
  # ax.set_xticks([1,3],["2023 December","2024 January"])
  ax.legend(loc="lower right",fontsize=20,ncol=2)
  ax.set_ylim(0,35)
  # ax.set_xlim(0,4)
  ax.grid(True,alpha=0.5)
  plt.savefig(plotFolder+"/triggerRateFRT2024Weekwise.pdf",transparent=False,bbox_inches='tight')
  plt.savefig(plotFolder+"/triggerRateFRT2024Weekwise.png",transparent=False,bbox_inches='tight')
