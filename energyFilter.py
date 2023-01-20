import matplotlib as mpl
mpl.use("Agg")
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
from matplotlib.ticker import AutoMinorLocator,MultipleLocator
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
from matplotlib import colors

from customColors import qualitative_colors

import numpy as np
import pickle

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

from weighting import GetWeight, ParticleType, PDGCode
from inclinedTriggerTools import *

from icecube.weighting.fluxes import GaisserH4a_IT


#run python reconstructionEfficiency.py --usePickle

basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/"



plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"


# colorsList = ['#9467bd', '#e377c2','#1f77b4','#2ca02c','#bcbd22','#ff7f0e','#8c564b','#7f7f7f','#17becf','#d62728',
#       '#4477AA', '#332288', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
#       '#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']
colorsList = ['#1f77b4','#ff7f0e','#2ca02c','#8c564b','#9467bd', '#e377c2','#bcbd22','#7f7f7f','#17becf','#d62728',
      '#4477AA', '#332288','#2ca02c', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
      '#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray',
      '#1f77b4','#ff7f0e','#2ca02c','#8c564b','#9467bd', '#e377c2','#bcbd22','#7f7f7f','#17becf','#d62728',
      '#4477AA', '#332288','#2ca02c', '#6699CC', '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77', '#661100', '#CC6677',
      '#AA4466','#882255','#AA4499'"tab:brown","lime",'teal',"navy","darkkhaki","olive","gray", 'sienna','slategray']

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)

# multiple_color = mpl.cm.tab20(range(20))
# multiple_color = np.concatenate((np.asarray([]),mpl.cm.winter(range(4)),mpl.cm.cool(range(4)),mpl.cm.Wistia(range(4)),mpl.cm.copper(range(4))))
multiple_color = ["olive","darkslategray","teal","darkturquoise","paleturquoise","midnightblue","royalblue","cornflowerblue","lightsteelblue","brown","indianred","lightcoral","rosybrown","darkgoldenrod","goldenrod","gold","khaki","purple","mediumorchid","violet","plum"]

triggerList = ["HLC6_5000","tank6_5000","tank6_4000","tank6_3000","tank6_2000",
  "tank7_5000","tank7_4000","tank7_3000","tank7_2000","tank8_5000","tank8_4000","tank8_3000","tank8_2000"]
triggerList7 = ["HLC6_5000","tank7_5000","tank7_4000","tank7_3000","tank7_2000"]
# triggerListSelect = ["HLC6_5000","tank6_3000","tank7_3000","tank8_3000"]
triggerListSelect = ["HLC6_5000","tank7_3000"]


trigWindow = 10**(-6) # in ns

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
# sin2ZenBins = [0.0,0.822]

energyBins = 10**np.linspace(5.0, 8.0, 31)
energyBinsShort = 10**np.linspace(5.0, 8.0, 7)
# energyBins = 10**np.linspace(5, 8.0, 7)
energyBinslgE = np.linspace(5.0,8.9,4000)
energyBinCenter = [5.1,6.1,7.1,8.1]

pickleFilesPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetPickle/"
pickleFiles = sorted(glob.glob(pickleFilesPath+"*GenDetFiltProcUniqueCleanVEMEvts.pkl"))


evtListReco = []
for ipickle in pickleFiles:
  with open(ipickle,"rb") as f:
    ievtListReco = pickle.load(f)
    evtListReco += ievtListReco

def plotHitsEnergy(eventList,hitType,triggerType,containment):
  if containment == True:
    eventList = containedEvents(eventList,410)
  if triggerType != None:
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  energy = [np.log10(ievt.energy*10**(9.0)) for ievt in eventList]
  hits = [getattr(ievt,hitType) for ievt in eventList]
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  # bins=[np.linspace(14,17,31),np.linspace(0,160,161)]
  if "SLC" in hitType:
    bins=[np.linspace(14,17,31),np.linspace(0,50,51)]
  else:
    bins=[np.linspace(14,17,31),np.linspace(0,130,131)]
  # ax.plot(energy,hits,".",c="orange",label=r"".format(hitType))
  counts,xedges,yedges,im=ax.hist2d(energy,hits,bins=bins,norm=colors.LogNorm(),label=r"".format(hitType))
  cbar = fig.colorbar(im, ax=ax)
  cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
  cbar.ax.set_ylabel(r"count", fontsize=22)
  ax.axvline(x=16.0,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,alpha=1.0)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  # ax.set_ylabel(r"trigger efficiency", fontsize=22)
  ax.set_ylabel(r"hits", fontsize=22)
  if triggerType != None:
    ax.set_title("{} {}".format(hitType,triggerType),fontsize=22)
  else:
    ax.set_title("{}".format(hitType),fontsize=22)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  # ax.text(0.78,0.1,s=r"trig:{0}".format(triggerType),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_xscale('log')
  # ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  # ax.set_ylim(0,40)
  # ax.set_xlim(14,17)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  ax.grid(True,alpha=0.5)
  # ax.legend(loc="upper left",fontsize=12)
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"0.98")
  # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  # ax.add_artist(l1)
  # ax.add_artist(l2)
  plt.savefig(plotFolder+"/hit{}EnergyCont{}Trig{}.pdf".format(hitType,containment,triggerType),transparent=False,bbox_inches='tight')
  plt.close()


# plotHitsEnergy(evtListReco,"nSLCTank","tank7_3000",True)
# plotHitsEnergy(evtListReco,"nSLCTank","tank7_3000",False)
# plotHitsEnergy(evtListReco,"nHLCTank","tank7_3000",True)
# plotHitsEnergy(evtListReco,"nHLCTank","tank7_3000",False)

# plotHitsEnergy(evtListReco,"nSLCTank",None,True)
# plotHitsEnergy(evtListReco,"nSLCTank",None,False)
# plotHitsEnergy(evtListReco,"nHLCTank",None,True)
# plotHitsEnergy(evtListReco,"nHLCTank",None,False)
# plotHitsEnergy(evtListReco,nHLC)

def plotHitsEnergyZenBin(eventList,hitType,triggerType,containment):
  if containment == True:
    eventList = containedEvents(eventList,410)
  if triggerType != None:
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
    energy = [np.log10(ievt.energy*10**(9.0)) for ievt in evtZenBin]
    hits = [getattr(ievt,hitType) for ievt in evtZenBin]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins=[np.linspace(14,17,31),np.linspace(0,160,161)]
    if "SLC" in hitType:
      bins=[np.linspace(14,17,31),np.linspace(0,50,51)]
    else:
      bins=[np.linspace(14,17,31),np.linspace(0,130,131)]
    # ax.plot(energy,hits,".",c="orange",label=r"".format(hitType))
    counts,xedges,yedges,im=ax.hist2d(energy,hits,bins=bins,norm=colors.LogNorm(),label=r"".format(hitType))
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
    cbar.ax.set_ylabel(r"count", fontsize=22)
    ax.axvline(x=16.0,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
    # ax.set_ylabel(r"trigger efficiency", fontsize=22)
    ax.set_ylabel(r"hits", fontsize=22)
    if triggerType != None:
      ax.set_title(r"{0} {1} {2:.1f}$^{{\circ}}$-{3:.1f}$^{{\circ}}$".format(hitType,triggerType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),fontsize=16)
    else:
      ax.set_title(r"{0} {1:.1f}$^{{\circ}}$-{2:.1f}$^{{\circ}}$".format(hitType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),fontsize=16)
    ax.grid(True,alpha=0.5)
    plt.savefig(plotFolder+"/hit{}EnergyCont{}Trig{}Zen{:.1f}_{:.1f}.pdf".format(
      hitType,containment,triggerType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),
    transparent=False,bbox_inches='tight')
    plt.close()


# plotHitsEnergyZenBin(evtListReco,"nSLCTank","tank7_3000",True)
# plotHitsEnergyZenBin(evtListReco,"nHLCTank","tank7_3000",True)






def plotQtotEnergy(eventList,QtotType,triggerType,containment):
  if containment == True:
    eventList = containedEvents(eventList,410)
  if triggerType != None:
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  energy = [np.log10(ievt.energy*10**(9.0)) for ievt in eventList]
  hits = [getattr(ievt,QtotType) for ievt in eventList]
  print("min,max",min(hits),max(hits),min(energy),max(energy))
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  # bins=[np.linspace(14,17,31),np.linspace(0,160,161)]
  # bins=[np.linspace(14,17,31),np.linspace(0,1200,201)]
  if "SLC" in QtotType:
    bins=[np.linspace(14,17,31),np.linspace(0,50,51)]
  else:
    bins=[np.linspace(14,17,31),np.linspace(0,1200,1201)]
  # ax.plot(energy,hits,".",c="orange",label=r"".format(hitType))
  # ax.hist2d(energy,hits,bins=bins,label=r"".format(QtotType))
  counts,xedges,yedges,im=ax.hist2d(energy,hits,bins=bins,norm=colors.LogNorm(),label=r"".format(QtotType))
  cbar = fig.colorbar(im, ax=ax)
  cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
  cbar.ax.set_ylabel(r"count", fontsize=22)
  ax.axvline(x=16.0,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,alpha=1.0)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  # ax.set_ylabel(r"trigger efficiency", fontsize=22)
  ax.set_ylabel(r"Qtot", fontsize=22)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  # ax.text(0.78,0.1,s=r"trig:{0}".format(triggerType),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_xscale('log')
  # ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  # ax.set_ylim(0,40)
  # ax.set_xlim(14,17)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  if triggerType != None:
    ax.set_title("{} {}".format(QtotType,triggerType),fontsize=22)
  else:
    ax.set_title("{}".format(QtotType),fontsize=22)
  ax.grid(True,alpha=0.5)
  # ax.legend(loc="upper left",fontsize=12)
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"0.98")
  # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  # ax.add_artist(l1)
  # ax.add_artist(l2)
  plt.savefig(plotFolder+"/Qtot{}EnergyCont{}Trig{}.pdf".format(QtotType,containment,triggerType),transparent=False,bbox_inches='tight')
  plt.close()


# plotQtotEnergy(evtListReco,hitType)
# plotQtotEnergy(evtListReco,"QtotSLCTank","tank7_3000",True)
# plotQtotEnergy(evtListReco,"QtotSLCTank","tank7_3000",False)
# plotQtotEnergy(evtListReco,"QtotHLCTank","tank7_3000",True)
# plotQtotEnergy(evtListReco,"QtotHLCTank","tank7_3000",False)

# plotQtotEnergy(evtListReco,"QtotSLCTank",None,True)
# plotQtotEnergy(evtListReco,"QtotSLCTank",None,False)
# plotQtotEnergy(evtListReco,"QtotHLCTank",None,True)
# plotQtotEnergy(evtListReco,"QtotHLCTank",None,False)


def plotQtotEnergyZenBin(eventList,QtotType,triggerType,containment):
  if containment == True:
    eventList = containedEvents(eventList,410)
  if triggerType != None:
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
    energy = [np.log10(ievt.energy*10**(9.0)) for ievt in evtZenBin]
    hits = [getattr(ievt,QtotType) for ievt in evtZenBin]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins=[np.linspace(14,17,31),np.linspace(0,160,161)]
    if "SLC" in QtotType:
      bins=[np.linspace(14,17,31),np.linspace(0,50,51)]
    else:
      bins=[np.linspace(14,17,31),np.linspace(0,200,201)]
    # ax.plot(energy,hits,".",c="orange",label=r"".format(QtotType))
    counts,xedges,yedges,im=ax.hist2d(energy,hits,bins=bins,norm=colors.LogNorm(),label=r"".format(QtotType))
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
    cbar.ax.set_ylabel(r"count", fontsize=22)
    ax.axvline(x=16.0,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
    # ax.set_ylabel(r"trigger efficiency", fontsize=22)
    ax.set_ylabel(r"Q$_{tot}$", fontsize=22)
    if triggerType != None:
      ax.set_title(r"{0} {1} {2:.1f}$^{{\circ}}$-{3:.1f}$^{{\circ}}$".format(QtotType,triggerType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),fontsize=16)
    else:
      ax.set_title(r"{0} {1:.1f}$^{{\circ}}$-{2:.1f}$^{{\circ}}$".format(QtotType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),fontsize=16)
    ax.grid(True,alpha=0.5)
    plt.savefig(plotFolder+"/hit{}EnergyCont{}Trig{}Zen{:.1f}_{:.1f}.pdf".format(
      QtotType,containment,triggerType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),
    transparent=False,bbox_inches='tight')
    plt.close()



# plotQtotEnergy(evtListReco,QtotType)
# plotQtotEnergyZenBin(evtListReco,"QtotSLCTank","tank7_3000",True)
plotQtotEnergyZenBin(evtListReco,"QtotHLCTank","tank7_3000",True)