#!/usr/bin/env python3

from icecube.icetray import I3Tray
from icecube import icetray,dataclasses,dataio

import re
import glob
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import datetime

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

runRateDict = {139170:24.2,139171:24.5,139172:24.42,139173:24.6,139174:25.33,139175:26.33,139176:26.55,
139177:26.33,139178:26.14,139179:25.68,139180:25.26,139181:24.88,139182:24.64,139183:24.52,139184:24.59,139185:24.75,139186:25.04,
139187:25.09,139188:25.33,139189:25.56,139190:25.89,139191:26.61,139192:27.66,139193:28.65,139194:29.05,139195:28.59,139196:27.02,
139197:25.57,139198:25.42,139199:24.87,139200:24.77,139201:25.2,139202:26.11,139203:26.93,139204:27.74,139205:28.33,139206:28.48,
139207:29.27,139208:29.67,139209:29.01,139210:28.3,139211:27.98,139212:27.6,139213:27.81,139214:28.29,139215:28.75,139216:29.07,
139217:29.02,139218:28.78,139219:28.42,139220:28.21,139221:28.13,139222:28.1,139223:27.98,139224:28.13}



fig = plt.figure(figsize=(8,5))
gs = gridspec.GridSpec(nrows=1,ncols=1)
ax = fig.add_subplot(gs[0])
ax.plot(runRateDict.keys(),runRateDict.values(),".",label="",alpha=1)
# ax.plot(runList,HLCrateList,".",label=r"ITSMT",alpha=1)
ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
secax = ax.secondary_xaxis('top')
secax.tick_params(axis='both',which='both', direction='in', labelsize=10)
ax.tick_params(axis='x', labelrotation=90)
secax.tick_params(axis='x', labelrotation=90)
# ax.set_xticks(runRateDict.keys())
# secax.set_xticks(runRateDict.keys())
# ax.xaxis.get_major_locator().set_params(integer=True)
# secax.set_xticklabels(dates)
# ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
# ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
ax.set_ylabel("rate [Hz]", fontsize=11)
# ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
ax.set_xlabel("Run", fontsize=11)
secax.set_xlabel("day", fontsize=11)
# ax.ticklabel_format(useOffset=False)
# ax.set_yscale("log")
# ax.set_xlim(-600,600)
# ax.set_ylim(-600,600)
# ax.legend()
plt.savefig(plotFolder+"/trigCountRunRateForbush.pdf",transparent=False,bbox_inches='tight')
plt.savefig(plotFolder+"/trigCountRunRateForbush.png",transparent=False,bbox_inches='tight')



from scipy import signal
fig = plt.figure(figsize=(8,5))
gs = gridspec.GridSpec(nrows=1,ncols=1)
ax = fig.add_subplot(gs[0])
values = list(runRateDict.values())
values = signal.detrend(values)
# mean = np.mean(values)
# values = [i-np.mean(values) for i in values]
ax.plot(runRateDict.keys(),values,".",label="",alpha=1)
# ax.plot(runList,HLCrateList,".",label=r"ITSMT",alpha=1)
ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
secax = ax.secondary_xaxis('top')
secax.tick_params(axis='both',which='both', direction='in', labelsize=10)
ax.tick_params(axis='x', labelrotation=90)
secax.tick_params(axis='x', labelrotation=90)
# ax.set_xticks(runRateDict.keys())
# secax.set_xticks(runRateDict.keys())
# ax.xaxis.get_major_locator().set_params(integer=True)
# secax.set_xticklabels(dates)
# ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
# ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
ax.set_ylabel("rate [Hz]", fontsize=11)
# ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
ax.set_xlabel("Run", fontsize=11)
secax.set_xlabel("day", fontsize=11)
# ax.ticklabel_format(useOffset=False)
# ax.set_yscale("log")
# ax.set_xlim(-600,600)
# ax.set_ylim(-600,600)
# ax.legend()
plt.savefig(plotFolder+"/trigCountRunRateForbushValues.pdf",transparent=False,bbox_inches='tight')
plt.savefig(plotFolder+"/trigCountRunRateForbushValues.png",transparent=False,bbox_inches='tight')