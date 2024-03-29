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

from scipy.optimize import curve_fit

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

def line(x,a,b):
  x = np.asarray(x)
  return a+b*x


def plotHitsEnergyZenBin(eventList,hitType,triggerType,containment,sin2ZenBins,xaxis,weighting,hitBin):
  if containment == True:
    eventList = containedEvents(eventList,410)
  if triggerType != None:
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
    if weighting == True:
      weights = [ievt.H4aWeight for ievt in evtZenBin]
    else:
      weights = [1 for ievt in evtZenBin]
    energy = [np.log10(ievt.energy*10**(9.0)) for ievt in evtZenBin]
    if "Diff" in hitType:
      hits = [getattr(ievt,"nSLCTank")-getattr(ievt,"nHLCTank") for ievt in evtZenBin]
    elif "Ratio" in hitType:
      hits = [getattr(ievt,"nSLCTank")/getattr(ievt,"nHLCTank") if getattr(ievt,"nHLCTank") > 0 else -100 for ievt in evtZenBin]
    elif "Sum" in hitType:
      hits = [getattr(ievt,"nSLCTank")+getattr(ievt,"nHLCTank") for ievt in evtZenBin]
    else:
      hits = [getattr(ievt,hitType) for ievt in evtZenBin]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins=[np.linspace(14,17,31),np.linspace(0,160,161)]
    eBin = np.linspace(14,17,31)
    # ax.plot(energy,hits,".",c="orange",label=r"".format(hitType))
    if xaxis == "energy":
      ylabel = r"hits" 
      xlabel = r"log10 (E [eV])"
      counts,xedges,yedges,im=ax.hist2d(energy,hits,bins=[eBin,hitBin],weights=weights,norm=colors.LogNorm(),label=r"".format(hitType))
    elif xaxis == "hits":
      xlabel = r"hits" 
      ylabel = r"log10 (E [eV])"
      counts,xedges,yedges,im=ax.hist2d(hits,energy,bins=[hitBin,eBin],weights=weights,norm=colors.LogNorm(),label=r"".format(hitType))
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
    cbar.ax.set_ylabel(r"count", fontsize=22)
    ################################################################
    # popt, pcov = curve_fit(line, hits, energy, p0=[14.0,1.5])
    # xfit = np.linspace(0,50,1000)
    # yfit = line(xfit,*popt)
    # ax.plot(xfit,yfit,ls="-",c="red",alpha=1.0)
    ################################################################
    ax.axhline(y=16.0,xmin=0,xmax=1,color="orange",ls="--",lw=2.5,alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    # ax.set_ylabel(r"trigger efficiency", fontsize=22)
    ax.set_xlabel(xlabel, fontsize=22)
    if triggerType != None:
      ax.set_title(r"{0} {1} {2:.1f}$^{{\circ}}$-{3:.1f}$^{{\circ}}$".format(hitType,triggerType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),fontsize=16)
    else:
      ax.set_title(r"{0} {1:.1f}$^{{\circ}}$-{2:.1f}$^{{\circ}}$".format(hitType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),fontsize=16)
    ax.grid(True,alpha=0.5)
    plt.savefig(plotFolder+"/hit{}EnergyCont{}Trig{}Zen{:.1f}_{:.1f}X{}.pdf".format(
      hitType,containment,triggerType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,xaxis),
    transparent=False,bbox_inches='tight')
    plt.close()


# plotHitsEnergyZenBin(evtListReco,"nSLCTank","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsEnergyZenBin(evtListReco,"nSLCTank","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsEnergyZenBin(evtListReco,"nSLCTank",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsEnergyZenBin(evtListReco,"nHLCTank","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsEnergyZenBin(evtListReco,"nHLCTank","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsEnergyZenBin(evtListReco,"nHLCTank",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))


plotHitsEnergyZenBin(evtListReco,"nSLCHLCSum","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin= np.linspace(0,150,151))
plotHitsEnergyZenBin(evtListReco,"nSLCHLCSum","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,150,151))
plotHitsEnergyZenBin(evtListReco,"nSLCHLCSum",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,150,151))


plotHitsEnergyZenBin(evtListReco,"nSLCHLCDiff","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin= np.linspace(-25,25,52))
plotHitsEnergyZenBin(evtListReco,"nSLCHLCDiff","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(-25,25,52))
plotHitsEnergyZenBin(evtListReco,"nSLCHLCDiff",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(-25,25,52))


plotHitsEnergyZenBin(evtListReco,"nSLCHLCRatio","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin= np.linspace(-1,15,151))
plotHitsEnergyZenBin(evtListReco,"nSLCHLCRatio","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(-1,15,151))
plotHitsEnergyZenBin(evtListReco,"nSLCHLCRatio",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(-1,15,151))

# plotHitsEnergyZenBin(evtListReco,"nSLCTank","tank7_3000",True,sin2ZenBins,xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsEnergyZenBin(evtListReco,"nSLCTank","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsEnergyZenBin(evtListReco,"nSLCTank",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsEnergyZenBin(evtListReco,"nHLCTank","tank7_3000",True,sin2ZenBins,xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsEnergyZenBin(evtListReco,"nHLCTank","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsEnergyZenBin(evtListReco,"nHLCTank",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))


def plotQtotEnergyZenBin(eventList,QtotType,triggerType,containment,sin2ZenBins,xaxis,weighting,hitBin):
  if containment == True:
    eventList = containedEvents(eventList,410)
  if triggerType != None:
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
    if weighting == True:
      weights = [ievt.H4aWeight for ievt in evtZenBin]
    else:
      weights = [1 for ievt in evtZenBin]
    energy = [np.log10(ievt.energy*10**(9.0)) for ievt in evtZenBin]
    if "Diff" in QtotType:
      hits = [getattr(ievt,"QtotSLCTank")-getattr(ievt,"QtotHLCTank") for ievt in evtZenBin]
    elif "Ratio" in QtotType:
      hits = [getattr(ievt,"QtotSLCTank")/getattr(ievt,"QtotHLCTank") if getattr(ievt,"QtotHLCTank")>0 else -100 for ievt in evtZenBin]
    elif "Sum" in QtotType:
      hits = [getattr(ievt,"QtotSLCTank")+getattr(ievt,"QtotHLCTank") for ievt in evtZenBin]
    else:
      hits = [getattr(ievt,QtotType) for ievt in evtZenBin]
    hits = [np.log10(ihit) if ihit > 0 else -10 for ihit in hits ]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins=[np.linspace(14,17,31),np.linspace(0,160,161)]
    eBin = np.linspace(14,17,31)
    # ax.plot(energy,hits,".",c="orange",label=r"".format(hitType))
    if xaxis == "energy":
      ylabel = r"log10 (Q$_{tot}$ [VEM])"
      xlabel = r"log10 (E [eV])"
      counts,xedges,yedges,im=ax.hist2d(energy,hits,bins=[eBin,hitBin],weights=weights,norm=colors.LogNorm(),label=r"".format(QtotType))
    elif xaxis == "hits":
      xlabel = r"log10 (Q$_{tot}$ [VEM])"
      ylabel = r"log10 (E [eV])"
      counts,xedges,yedges,im=ax.hist2d(hits,energy,bins=[hitBin,eBin],weights=weights,norm=colors.LogNorm(),label=r"".format(QtotType))
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
    cbar.ax.set_ylabel(r"count", fontsize=22)
    ax.axhline(y=16.0,xmin=0,xmax=1,color="orange",ls="--",lw=2.5,alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    # ax.set_ylabel(r"trigger efficiency", fontsize=22)
    ax.set_xlabel(xlabel, fontsize=22)
    if triggerType != None:
      ax.set_title(r"{0} {1} {2:.1f}$^{{\circ}}$-{3:.1f}$^{{\circ}}$".format(QtotType,triggerType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),fontsize=16)
    else:
      ax.set_title(r"{0} {1:.1f}$^{{\circ}}$-{2:.1f}$^{{\circ}}$".format(QtotType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),fontsize=16)
    ax.grid(True,alpha=0.5)
    plt.savefig(plotFolder+"/hit{}EnergyCont{}Trig{}Zen{:.1f}_{:.1f}X{}.pdf".format(
      QtotType,containment,triggerType,np.arcsin(np.sqrt(sin2ZenBins[nbin])
        )*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,xaxis),
    transparent=False,bbox_inches='tight')
    plt.close()



# plotQtotEnergyZenBin(evtListReco,"QtotSLCTank","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin= np.linspace(0,200,201))
# plotQtotEnergyZenBin(evtListReco,"QtotSLCTank","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,200,201))
# plotQtotEnergyZenBin(evtListReco,"QtotSLCTank",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,200,201))
# plotQtotEnergyZenBin(evtListReco,"QtotHLCTank","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin= np.linspace(0,200,201))
# plotQtotEnergyZenBin(evtListReco,"QtotHLCTank","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,200,201))
# plotQtotEnergyZenBin(evtListReco,"QtotHLCTank",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,200,201))


plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCSum","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin=np.linspace(-1,4,51))
plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCSum","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-1,4,51))
plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCSum",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-1,4,51))

plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCDiff","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin=np.linspace(-5,5,51))
plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCDiff","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-5,5,51))
plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCDiff",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-5,5,51))


plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCRatio","tank7_3000",True,sin2ZenBins,xaxis="hits",weighting=True,hitBin=np.linspace(-5,5,51))
plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCRatio","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-5,5,51))
plotQtotEnergyZenBin(evtListReco,"QtotSLCHLCRatio",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-5,5,51))


# plotQtotEnergyZenBin(evtListReco,"QtotSLCTank","tank7_3000",True,sin2ZenBins,xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotEnergyZenBin(evtListReco,"QtotSLCTank","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotEnergyZenBin(evtListReco,"QtotSLCTank",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotEnergyZenBin(evtListReco,"QtotHLCTank","tank7_3000",True,sin2ZenBins,xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotEnergyZenBin(evtListReco,"QtotHLCTank","tank7_3000",True,[sin2ZenBins[i] for i in [0,-1]],xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotEnergyZenBin(evtListReco,"QtotHLCTank",None,True,[sin2ZenBins[i] for i in [0,-1]],xaxis="energy",weighting=True,hitBin= np.linspace(0,300,301))

##########################################################################################
##########################################################################################
################as a function of zenith angle ############################################
def plotHitsZenithEnergyBin(eventList,hitType,triggerType,containment,energyBins,xaxis,weighting,hitBin):
  if containment == True:
    eventList = containedEvents(eventList,410)
  if triggerType != None:
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  for nbin, binStart in enumerate(energyBins[:-1]):
    lowEdge = energyBins[nbin]
    highEdge = energyBins[nbin+1]
    evtEBin = [ievt for ievt in eventList if lowEdge <= ievt.energy < highEdge]
    if weighting == True:
      weights = [ievt.H4aWeight for ievt in evtEBin]
    else:
      weights = [1 for ievt in evtEBin]
    # zenith = [ievt.zenith*180.0/np.pi for ievt in evtEBin]
    zenith = [np.sin(ievt.zenith)**2 for ievt in evtEBin]
    if "Diff" in hitType:
      hits = [getattr(ievt,"nSLCTank")-getattr(ievt,"nHLCTank") for ievt in evtEBin]
    elif "Ratio" in hitType:
      hits = [getattr(ievt,"nSLCTank")/getattr(ievt,"nHLCTank") if getattr(ievt,"nHLCTank") > 0 else -100 for ievt in evtEBin]
    elif "Sum" in hitType:
      hits = [getattr(ievt,"nSLCTank")+getattr(ievt,"nHLCTank") for ievt in evtEBin]
    else:
      hits = [getattr(ievt,hitType) for ievt in evtEBin]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins=[np.linspace(14,17,31),np.linspace(0,160,161)]
    # zenBin = np.linspace(0,66,67)
    zenBin = sin2ZenBins
    # ax.plot(energy,hits,".",c="orange",label=r"".format(hitType))
    if xaxis == "zenith":
      # xlabel = r"$\theta^{\circ}$"
      xlabel = r"sin $^2\,\theta$"
      ylabel = r"hits"
      counts,xedges,yedges,im=ax.hist2d(zenith,hits,bins=[zenBin,hitBin],weights=weights,norm=colors.LogNorm(),label=r"".format(hitType))
    elif xaxis == "hits":
      # ylabel = r"$\theta^{\circ}$"
      ylabel = r"sin $^2\,\theta$"
      xlabel = r"hits"
      counts,xedges,yedges,im=ax.hist2d(hits,zenith,bins=[hitBin,zenBin],weights=weights,norm=colors.LogNorm(),label=r"".format(hitType))
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
    cbar.ax.set_ylabel(r"count", fontsize=22)
    # ax.axhline(y=16.0,xmin=0,xmax=1,color="orange",ls="--",lw=2.5,alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    # ax.set_ylabel(r"trigger efficiency", fontsize=22)
    ax.set_xlabel(xlabel, fontsize=22)
    if triggerType != None:
      ax.set_title(r"{0} {1} lgE [eV] = {2:.1f}-{3:.1f}".format(hitType,triggerType,np.log10(energyBins[nbin])+9
        ,np.log10(energyBins[nbin+1])+9),fontsize=16)
    else:
      ax.set_title(r"{0} lgE [eV] = {1:.1f}-{2:.1f}".format(hitType,np.log10(energyBins[nbin])+9
        ,np.log10(energyBins[nbin+1])+9),fontsize=16)
    ax.grid(True,alpha=0.5)
    plt.savefig(plotFolder+"/hit{}ZenithCont{}Trig{}Energy{:.1f}_{:.1f}X{}.pdf".format(
      hitType,containment,triggerType,np.log10(energyBins[nbin])+9
        ,np.log10(energyBins[nbin+1])+9,xaxis),transparent=False,bbox_inches='tight')
    plt.close()


# plotHitsZenithEnergyBin(evtListReco,"nSLCTank","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsZenithEnergyBin(evtListReco,"nSLCTank","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsZenithEnergyBin(evtListReco,"nSLCTank",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsZenithEnergyBin(evtListReco,"nHLCTank","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsZenithEnergyBin(evtListReco,"nHLCTank","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotHitsZenithEnergyBin(evtListReco,"nHLCTank",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))



plotHitsZenithEnergyBin(evtListReco,"nSLCHLCSum","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin= np.linspace(0,150,151))
plotHitsZenithEnergyBin(evtListReco,"nSLCHLCSum","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,150,151))
plotHitsZenithEnergyBin(evtListReco,"nSLCHLCSum",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,150,151))


plotHitsZenithEnergyBin(evtListReco,"nSLCHLCDiff","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin= np.linspace(-70,25,96))
plotHitsZenithEnergyBin(evtListReco,"nSLCHLCDiff","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(-70,25,96))
plotHitsZenithEnergyBin(evtListReco,"nSLCHLCDiff",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(-70,25,96))


plotHitsZenithEnergyBin(evtListReco,"nSLCHLCRatio","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin= np.linspace(-1,14,151))
plotHitsZenithEnergyBin(evtListReco,"nSLCHLCRatio","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(-1,14,151))
plotHitsZenithEnergyBin(evtListReco,"nSLCHLCRatio",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(-1,14,151))


# plotHitsZenithEnergyBin(evtListReco,"nSLCTank","tank7_3000",True,energyBinsShort,xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsZenithEnergyBin(evtListReco,"nSLCTank","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsZenithEnergyBin(evtListReco,"nSLCTank",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsZenithEnergyBin(evtListReco,"nHLCTank","tank7_3000",True,energyBinsShort,xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsZenithEnergyBin(evtListReco,"nHLCTank","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotHitsZenithEnergyBin(evtListReco,"nHLCTank",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))


def plotQtotZenithEnergyBin(eventList,QtotType,triggerType,containment,energyBins,xaxis,weighting,hitBin):
  if containment == True:
    eventList = containedEvents(eventList,410)
  if triggerType != None:
    eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  for nbin, binStart in enumerate(energyBins[:-1]):
    lowEdge = energyBins[nbin]
    highEdge = energyBins[nbin+1]
    evtEBin = [ievt for ievt in eventList if lowEdge <= ievt.energy < highEdge]
    if weighting == True:
      weights = [ievt.H4aWeight for ievt in evtEBin]
    else:
      weights = [1 for ievt in evtEBin]
    # zenith = [ievt.zenith*180.0/np.pi for ievt in evtEBin]
    zenith = [np.sin(ievt.zenith)**2 for ievt in evtEBin]
    if "Diff" in QtotType:
      hits = [getattr(ievt,"QtotSLCTank")-getattr(ievt,"QtotHLCTank") for ievt in evtEBin]
    elif "Ratio" in QtotType:
      hits = [getattr(ievt,"QtotSLCTank")/getattr(ievt,"QtotHLCTank") if getattr(ievt,"QtotHLCTank") > 0 else -100 for ievt in evtEBin]
    elif "Sum" in QtotType:
      hits = [getattr(ievt,"QtotSLCTank")+getattr(ievt,"QtotHLCTank") for ievt in evtEBin]
    else:
      hits = [getattr(ievt,QtotType) for ievt in evtEBin]
    hits = [np.log10(ihit) if ihit > 0 else -10 for ihit in hits ]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins=[np.linspace(14,17,31),np.linspace(0,160,161)]
    # zenBin = np.linspace(np.sin(0)**2,np.sin(65*np.pi/180.0)**2,)
    zenBin = sin2ZenBins
    # ax.plot(energy,hits,".",c="orange",label=r"".format(QtotType))
    if xaxis == "zenith":
      # xlabel = r"$\theta^{\circ}$"
      xlabel = r"sin $^2\,\theta$"
      ylabel = r"log10(Q$_{tot}$ [VEM])"
      yscale = "log"
      xscale = "linear"
      counts,xedges,yedges,im=ax.hist2d(zenith,hits,bins=[zenBin,hitBin],weights=weights,norm=colors.LogNorm(),label=r"".format(QtotType))
    elif xaxis == "hits":
      # ylabel = r"$\theta^{\circ}$"
      ylabel = r"sin $^2\,\theta$"
      xlabel = r"log10(Q$_{tot}$ [VEM])"
      xscale = "log"
      yscale = "linear"
      counts,xedges,yedges,im=ax.hist2d(hits,zenith,bins=[hitBin,zenBin],weights=weights,norm=colors.LogNorm(),label=r"".format(QtotType))
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    cbar.set_ticks(mtick.LogLocator(), update_ticks=True)
    cbar.ax.set_ylabel(r"count", fontsize=22)
    # ax.axhline(y=16.0,xmin=0,xmax=1,color="orange",ls="--",lw=2.5,alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    # ax.set_ylabel(r"trigger efficiency", fontsize=22)
    ax.set_xlabel(xlabel, fontsize=22)
    if triggerType != None:
      ax.set_title(r"{0} {1} lgE [eV] = {2:.1f}-{3:.1f}".format(QtotType,triggerType,np.log10(energyBins[nbin])+9
        ,np.log10(energyBins[nbin+1])+9),fontsize=16)
    else:
      ax.set_title(r"{0} lgE [eV] = {1:.1f}-{2:.1f}".format(QtotType,np.log10(energyBins[nbin])+9
        ,np.log10(energyBins[nbin+1])+9),fontsize=16)
    ax.grid(True,alpha=0.5)
    plt.savefig(plotFolder+"/hit{}ZenithCont{}Trig{}E{:.1f}_{:.1f}X{}.pdf".format(
      QtotType,containment,triggerType,np.log10(energyBins[nbin])+9
        ,np.log10(energyBins[nbin+1])+9,xaxis),transparent=False,bbox_inches='tight')
    plt.close()



# plotQtotEnergy(evtListReco,QtotType)
# plotQtotZenithEnergyBin(evtListReco,"QtotSLCTank","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotQtotZenithEnergyBin(evtListReco,"QtotSLCTank","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotQtotZenithEnergyBin(evtListReco,"QtotSLCTank",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,50,51))
# plotQtotZenithEnergyBin(evtListReco,"QtotHLCTank","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotZenithEnergyBin(evtListReco,"QtotHLCTank","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotZenithEnergyBin(evtListReco,"QtotHLCTank",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin= np.linspace(0,300,301))


plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCSum","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin=np.linspace(-1,4,51))
plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCSum","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-1,4,51))
plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCSum",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-1,4,51))


plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCDiff","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin=np.linspace(-2.5,2.1,51))
plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCDiff","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-2.5,2.1,51))
plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCDiff",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-2.5,2.1,51))

plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCRatio","tank7_3000",True,energyBinsShort,xaxis="hits",weighting=True,hitBin=np.linspace(- 3,2,51))
plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCRatio","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-3,2,51))
plotQtotZenithEnergyBin(evtListReco,"QtotSLCHLCRatio",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="hits",weighting=True,hitBin=np.linspace(-3,2,51))


# plotQtotZenithEnergyBin(evtListReco,"QtotSLCTank","tank7_3000",True,energyBinsShort,xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotZenithEnergyBin(evtListReco,"QtotSLCTank","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotZenithEnergyBin(evtListReco,"QtotSLCTank",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotZenithEnergyBin(evtListReco,"QtotHLCTank","tank7_3000",True,energyBinsShort,xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotZenithEnergyBin(evtListReco,"QtotHLCTank","tank7_3000",True,[energyBinsShort[i] for i in [0,-1]],xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))
# plotQtotZenithEnergyBin(evtListReco,"QtotHLCTank",None,True,[energyBinsShort[i] for i in [0,-1]],xaxis="zenith",weighting=True,hitBin= np.linspace(0,300,301))