#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import datetime

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

with open("/home/enpaudel/icecube/triggerStudy/triggerRateList.txt","r") as f:
  run = []
  sub_run = []
  HG7rate = []
  HLCrate = []
  for line in f:
    splits = line.split(" ")
    if len(splits) > 3:
      run.append(int(splits[0]))
      sub_run.append(int(splits[1]))
      HG7rate.append(int(splits[2]))
      HLCrate.append(int(splits[3]))
runList = sorted(list(set(run)),reverse=True)
print("runList",runList)
rateDict = {}
rateDictHLC = {}
for irun in runList:
  rateDict[irun] = []
  rateDictHLC[irun] = []
for irun,irate in zip(run,HG7rate):
  rateDict[irun].append(irate)

for irun,irate in zip(run,HLCrate):
  rateDictHLC[irun].append(irate)



durationDict = {}
durationList = [7*60*60+2*60+58,8*60*60+60+3,8*60*60+60+8,8*60*60+60+11,
8*60*60+60+10,8*60*60+60+1,6*60*60+51*60+28,8*60*60+60+9,8*60*60+60+10,
7*60*60+59*60+59,8*60*60+60+10,8*60*60+60+11,8*60*60+60+7,
7*60*60+59*60+56,8*60*60+60+12,8*60*60+60+12,8*60*60+60+11,8*60*60+60+1]
for irun,iduration in zip(runList,durationList):
  durationDict[irun] = iduration
print("durationDict",durationDict)
HLCrateList = []
HG7rateList = []
# for irun in rateDict:
#   HLCrateList.append(sum(rateDictHLC[irun])/durationDict[irun])
#   HG7rateList.append(sum(rateDict[irun])/durationDict[irun])

# fig = plt.figure(figsize=(8,5))
# gs = gridspec.GridSpec(nrows=1,ncols=1)
# ax = fig.add_subplot(gs[0])
# ax.plot(runList,HG7rateList,".",label=r"IT7HG",alpha=1)
# ax.plot(runList,HLCrateList,".",label=r"ITSMT",alpha=1)
# ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
# ax.tick_params(axis='x', labelrotation=90)
# ax.set_xticks(runList)
# # ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
# # ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
# ax.set_ylabel("rate [Hz]", fontsize=22)
# # ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
# ax.set_xlabel("Run", fontsize=11)
# # ax.set_yscale("log")
# # ax.set_xlim(-600,600)
# # ax.set_ylim(-600,600)
# ax.legend()
# plt.savefig(plotFolder+"/trigCount.pdf",transparent=False,bbox_inches='tight')
# plt.savefig(plotFolder+"/hitsCount.png",transparent=False,bbox_inches='tight')

####################################################################################
####################################################################################
import pandas as pd
import numpy as np
GRL = "/data/exp/IceCube/2023/filtered/level2/"
GRL2023 = GRL + "IC86_2023_GoodRunInfo.txt"
GRL2023_df = pd.read_csv(GRL2023,sep="\\s+",header=0,escapechar='#')
GRL_IT_List = np.asarray(GRL2023_df.iloc[:,0][GRL2023_df.iloc[:,2]==1])
liveTime = np.asarray(GRL2023_df.iloc[:,3][GRL2023_df.iloc[:,2]==1])
runList = [int(irun) for irun in GRL_IT_List]
# runList = np.linspace(138935,138952,18)
HG7rateList = []
HLCrateList = []
for run,idur in zip(runList,liveTime):  
  # with open("/home/enpaudel/icecube/triggerStudy/triggerRateListCombinedShort{}.txt".format(int(run)),"r") as f:
  with open("/home/enpaudel/icecube/triggerStudy/rateFiles/triggerRate{}.txt".format(int(run)),"r") as f:
    iHG7rate = []
    iHLCrate = []
    for line in f:
      splits = line.split(" ")
      if len(splits) > 2:
        iHG7rate.append(int(splits[2]))
        iHLCrate.append(int(splits[3]))
  # HG7rateList.append(sum(iHG7rate)/durationDict[run])
  # HLCrateList.append(sum(iHLCrate)/durationDict[run])
  HG7rateList.append(sum(iHG7rate)/idur)
  HLCrateList.append(sum(iHLCrate)/idur)
  # print("run {:d} HG7 {:.2f} HLC6 {:.2f}".format(run,sum(iHG7rate)/idur,sum(iHLCrate)/idur))
dates = [np.datetime64("2023-11-28"),np.datetime64("2023-12-01"),np.datetime64("2023-12-04"),
np.datetime64("2023-12-07"),np.datetime64("2023-12-10"),np.datetime64("2023-12-12"),
np.datetime64("2023-12-16"),np.datetime64("2023-12-19"),np.datetime64("2023-12-21"),
np.datetime64("2023-12-23"),np.datetime64("2023-12-26"),np.datetime64("2023-12-28"),np.datetime64("2023-12-30"),
np.datetime64("2024-01-01"),np.datetime64("2024-01-05"),np.datetime64("2024-01-07"),
np.datetime64("2024-01-10"),np.datetime64("2024-01-13"),np.datetime64("2024-01-17"),np.datetime64("2024-01-20")]
# dates = [np.datetime64(iday) for iday in dates]

fig = plt.figure(figsize=(8,5))
gs = gridspec.GridSpec(nrows=1,ncols=1)
ax = fig.add_subplot(gs[0])
ax.plot(runList,HG7rateList,".",label=r"IT7HG",alpha=1)
ax.plot(runList,HLCrateList,".",label=r"ITSMT",alpha=1)
ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
secax = ax.secondary_xaxis('top')
secax.tick_params(axis='both',which='both', direction='in', labelsize=10)
ax.tick_params(axis='x', labelrotation=90)
secax.tick_params(axis='x', labelrotation=90)
ax.set_xticks(runList[::10])
secax.set_xticks(runList[::10])
secax.set_xticklabels(dates)
# ax.set_xlabel(r"$Q_{tot}$ [m]", fontsize=22)
# ax.set_ylabel(r"$\theta$$^{\circ}$", fontsize=22)
ax.set_ylabel("rate [Hz]", fontsize=11)
# ax.set_xlabel("Q$_{tot} [VEM]$", fontsize=22)
ax.set_xlabel("Run", fontsize=11)
secax.set_xlabel("day", fontsize=11)
# ax.set_yscale("log")
# ax.set_xlim(-600,600)
# ax.set_ylim(-600,600)
ax.legend()
plt.savefig(plotFolder+"/trigCount.pdf",transparent=False,bbox_inches='tight')
plt.savefig(plotFolder+"/trigCount.png",transparent=False,bbox_inches='tight')



# rates
# run 138615 HG7 25.42 HLC6 6.52
# run 138616 HG7 25.36 HLC6 6.53
# run 138617 HG7 25.65 HLC6 6.60
# run 138618 HG7 26.12 HLC6 6.67
# run 138619 HG7 26.20 HLC6 6.77
# run 138620 HG7 26.25 HLC6 6.68
# run 138621 HG7 26.36 HLC6 6.72
# run 138622 HG7 26.35 HLC6 6.73
# run 138623 HG7 26.47 HLC6 6.74
# run 138624 HG7 26.39 HLC6 6.73
# run 138626 HG7 26.50 HLC6 6.77
# run 138627 HG7 26.48 HLC6 6.75
# run 138628 HG7 25.92 HLC6 6.66
# run 138630 HG7 24.99 HLC6 6.54
# run 138631 HG7 24.88 HLC6 6.52
# run 138632 HG7 25.02 HLC6 6.48
# run 138633 HG7 24.85 HLC6 6.46
# run 138634 HG7 24.93 HLC6 6.44
# run 138635 HG7 25.18 HLC6 6.48
# run 138636 HG7 25.59 HLC6 6.57
# run 138637 HG7 26.48 HLC6 6.73
# run 138638 HG7 27.65 HLC6 7.02
# run 138647 HG7 28.33 HLC6 7.18
# run 138648 HG7 27.84 HLC6 7.05
# run 138649 HG7 26.47 HLC6 6.72
# run 138650 HG7 25.44 HLC6 6.54
# run 138651 HG7 25.46 HLC6 6.54
# run 138652 HG7 26.18 HLC6 6.70
# run 138653 HG7 27.10 HLC6 6.93
# run 138654 HG7 27.33 HLC6 6.97
# run 138655 HG7 27.33 HLC6 6.97
# run 138656 HG7 27.07 HLC6 6.91
# run 138657 HG7 7.72 HLC6 3.88
# run 138658 HG7 4.91 HLC6 2.12
# run 138660 HG7 26.80 HLC6 6.86
# run 138661 HG7 26.71 HLC6 6.83
# run 138662 HG7 26.96 HLC6 6.92
# run 138663 HG7 27.29 HLC6 6.95
# run 138664 HG7 27.62 HLC6 7.02
# run 138665 HG7 27.91 HLC6 7.06
# run 138666 HG7 28.44 HLC6 7.20
# run 138667 HG7 29.00 HLC6 7.28
# run 138668 HG7 28.76 HLC6 7.22
# run 138669 HG7 27.99 HLC6 7.09
# run 138670 HG7 27.55 HLC6 6.97
# run 138671 HG7 27.55 HLC6 6.95
# run 138674 HG7 27.69 HLC6 6.89
# run 138675 HG7 27.43 HLC6 7.05
# run 138677 HG7 27.43 HLC6 6.94
# run 138678 HG7 27.50 HLC6 6.98
# run 138680 HG7 27.58 HLC6 6.98
# run 138681 HG7 27.63 HLC6 7.03
# run 138682 HG7 27.23 HLC6 6.92
# run 138683 HG7 26.74 HLC6 6.82
# run 138684 HG7 26.62 HLC6 6.83
# run 138685 HG7 26.77 HLC6 6.89
# run 138686 HG7 26.82 HLC6 6.90
# run 138687 HG7 26.82 HLC6 6.90
# run 138688 HG7 26.44 HLC6 6.83
# run 138689 HG7 26.13 HLC6 6.77
# run 138690 HG7 25.97 HLC6 6.75
# run 138691 HG7 25.75 HLC6 6.70
# run 138692 HG7 26.23 HLC6 6.76
# run 138693 HG7 26.80 HLC6 6.86
# run 138694 HG7 27.35 HLC6 7.01
# run 138695 HG7 27.66 HLC6 7.13
# run 138696 HG7 33.53 HLC6 8.70
# run 138697 HG7 27.69 HLC6 7.13
# run 138698 HG7 27.25 HLC6 7.07
# run 138699 HG7 27.10 HLC6 7.03
# run 138700 HG7 27.31 HLC6 7.07
# run 138701 HG7 27.36 HLC6 7.11
# run 138702 HG7 27.36 HLC6 7.06
# run 138703 HG7 27.23 HLC6 7.07
# run 138704 HG7 27.07 HLC6 7.03
# run 138706 HG7 26.92 HLC6 7.00
# run 138707 HG7 26.92 HLC6 7.06
# run 138708 HG7 26.77 HLC6 6.97
# run 138709 HG7 26.83 HLC6 6.98
# run 138710 HG7 26.59 HLC6 6.85
# run 138711 HG7 26.61 HLC6 6.93
# run 138729 HG7 25.90 HLC6 6.68
# run 138734 HG7 4.95 HLC6 2.12
# run 138735 HG7 26.18 HLC6 6.91
# run 138739 HG7 26.47 HLC6 6.85
# run 138740 HG7 26.43 HLC6 6.89
# run 138741 HG7 26.34 HLC6 6.86
# run 138742 HG7 26.12 HLC6 6.80
# run 138743 HG7 25.98 HLC6 6.81
# run 138751 HG7 26.04 HLC6 6.78
# run 138752 HG7 26.21 HLC6 6.87
# run 138753 HG7 26.25 HLC6 6.87
# run 138754 HG7 26.71 HLC6 6.92
# run 138755 HG7 27.35 HLC6 7.03
# run 138756 HG7 28.12 HLC6 7.17
# run 138757 HG7 28.90 HLC6 7.36
# run 138758 HG7 29.32 HLC6 7.44
# run 138759 HG7 29.44 HLC6 7.45
# run 138760 HG7 29.33 HLC6 7.39
# run 138761 HG7 29.03 HLC6 7.35
# run 138762 HG7 28.94 HLC6 7.34
# run 138763 HG7 28.65 HLC6 7.28
# run 138764 HG7 12.58 HLC6 4.58
# run 138765 HG7 28.52 HLC6 7.27
# run 138766 HG7 28.35 HLC6 7.19
# run 138767 HG7 28.32 HLC6 7.19
# run 138768 HG7 2.27 HLC6 1.59
# run 138769 HG7 28.47 HLC6 7.24
# run 138770 HG7 28.50 HLC6 7.27
# run 138771 HG7 28.56 HLC6 7.27
# run 138772 HG7 27.88 HLC6 7.18
# run 138773 HG7 27.40 HLC6 7.06
# run 138774 HG7 26.74 HLC6 6.89
# run 138782 HG7 26.38 HLC6 6.94
# run 138783 HG7 26.26 HLC6 6.69
# run 138784 HG7 26.71 HLC6 6.91
# run 138786 HG7 27.16 HLC6 7.02
# run 138787 HG7 26.37 HLC6 6.85
# run 138788 HG7 26.21 HLC6 6.83
# run 138793 HG7 26.02 HLC6 6.74
# run 138794 HG7 26.10 HLC6 6.76
# run 138801 HG7 25.89 HLC6 6.80
# run 138802 HG7 25.89 HLC6 6.77
# run 138803 HG7 25.89 HLC6 6.78
# run 138804 HG7 25.78 HLC6 6.76
# run 138805 HG7 25.70 HLC6 6.76
# run 138806 HG7 25.79 HLC6 6.76
# run 138807 HG7 25.87 HLC6 6.74
# run 138808 HG7 26.18 HLC6 6.86
# run 138809 HG7 26.45 HLC6 6.90
# run 138810 HG7 26.77 HLC6 6.93
# run 138811 HG7 26.77 HLC6 6.92
# run 138812 HG7 26.61 HLC6 6.93
# run 138813 HG7 26.59 HLC6 6.89
# run 138814 HG7 26.39 HLC6 6.83
# run 138815 HG7 26.24 HLC6 6.83
# run 138816 HG7 26.19 HLC6 6.78
# run 138817 HG7 25.88 HLC6 6.74
# run 138818 HG7 25.94 HLC6 6.78
# run 138819 HG7 26.09 HLC6 6.80
# run 138820 HG7 26.25 HLC6 6.84
# run 138821 HG7 26.31 HLC6 6.86
# run 138822 HG7 26.45 HLC6 6.87
# run 138823 HG7 7.54 HLC6 3.83
# run 138824 HG7 4.85 HLC6 2.10
# run 138826 HG7 26.37 HLC6 6.72
# run 138827 HG7 26.32 HLC6 6.85
# run 138828 HG7 26.24 HLC6 6.80
# run 138829 HG7 26.09 HLC6 6.79
# run 138830 HG7 25.94 HLC6 6.76
# run 138831 HG7 26.01 HLC6 6.78
# run 138832 HG7 25.97 HLC6 6.74
# run 138833 HG7 26.04 HLC6 6.74
# run 138834 HG7 26.19 HLC6 6.80
# run 138835 HG7 26.38 HLC6 6.82
# run 138836 HG7 26.39 HLC6 6.84
# run 138837 HG7 26.27 HLC6 6.81
# run 138847 HG7 25.77 HLC6 6.68
# run 138848 HG7 25.70 HLC6 6.77
# run 138849 HG7 25.62 HLC6 6.70
# run 138850 HG7 25.72 HLC6 6.70
# run 138851 HG7 25.83 HLC6 6.71
# run 138852 HG7 26.11 HLC6 6.76
# run 138853 HG7 26.52 HLC6 6.86
# run 138854 HG7 26.84 HLC6 6.93
# run 138855 HG7 26.70 HLC6 6.90
# run 138856 HG7 26.41 HLC6 6.85
# run 138857 HG7 25.86 HLC6 6.76
# run 138858 HG7 25.29 HLC6 6.65
# run 138859 HG7 25.08 HLC6 6.56
# run 138860 HG7 25.24 HLC6 6.63
# run 138861 HG7 25.45 HLC6 6.68
# run 138862 HG7 25.59 HLC6 6.68
# run 138863 HG7 25.49 HLC6 6.66
# run 138864 HG7 25.38 HLC6 6.66
# run 138865 HG7 25.30 HLC6 6.62
# run 138866 HG7 25.01 HLC6 6.58
# run 138867 HG7 24.78 HLC6 6.55
# run 138868 HG7 24.51 HLC6 6.46
# run 138869 HG7 24.30 HLC6 6.42
# run 138870 HG7 24.33 HLC6 6.42
# run 138871 HG7 24.53 HLC6 6.46
# run 138872 HG7 24.64 HLC6 6.50
# run 138873 HG7 24.57 HLC6 6.46
# run 138874 HG7 24.61 HLC6 6.49
# run 138875 HG7 24.38 HLC6 6.41
# run 138876 HG7 24.16 HLC6 6.35
# run 138877 HG7 24.03 HLC6 6.33
# run 138878 HG7 24.12 HLC6 6.34
# run 138879 HG7 24.45 HLC6 6.45
# run 138880 HG7 24.86 HLC6 6.52
# run 138881 HG7 25.12 HLC6 6.55
# run 138882 HG7 25.16 HLC6 6.58
# run 138883 HG7 25.30 HLC6 6.65
