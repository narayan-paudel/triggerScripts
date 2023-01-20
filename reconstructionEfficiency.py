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

print("reading the evtObj")
evtList = []
for ipickle in pickleFiles:
  with open(ipickle,"rb") as f:
    ievtList = pickle.load(f)
    evtList += ievtList
idList = [ievt.eventID for ievt in evtList]
print("evtlsit",len(idList),min(idList),max(idList))

def plotRecoEfficiency(evtList,energyBins,triggerType,containment,recoThreshold,denoType):
  '''
  plots reco efficiency in different zenith bins
  '''
  print("plotting reco efficiency for ",triggerType)
  if containment == True:
    # evtList = containedEvents(evtList,640)
    evtList = containedEvents(evtList,410)
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # colorIter = iter(colorsList)
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
    energyList = []
    efficiencyList = []
    recoEfficiencyList = []
    ncolor = colorsCustom2[nbin]
    for ebin, ebinStart in enumerate(energyBins[:-1]):
      lowEdge_E = energyBins[ebin]
      highEdge_E = energyBins[ebin+1]
      evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
      # totalEvts = len(evtEBin)
      weights = [ievt.H4aWeight for ievt in evtEBin]
      totalEvts = len(evtEBin)
      # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
      # triggerList = [ievt for ievt in evtEBin if abs(getattr(ievt,triggerType)-1)<0.01 and abs(getattr(ievt,"HLC6_5000")-0)<0.01]
      triggerList = [ievt for ievt in evtEBin if abs(getattr(ievt,triggerType)-1)<0.01]
      # recoList = [ievt for ievt in triggerList if abs(ievt.zenithReco) <= np.pi]
      recoList = [ievt for ievt in triggerList if abs(ievt.zenithReco-ievt.zenith)*180.0/np.pi <= recoThreshold]
      trigEff = triggerEfficiency(len(triggerList),totalEvts)
      if denoType == "trig":
        recoEff = recoEfficiency(len(recoList),len(triggerList))
      elif denoType == "total":
        recoEff = recoEfficiency(len(recoList),totalEvts)
      efficiencyList.append(trigEff)
      recoEfficiencyList.append(recoEff)
      energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))   
    # ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
    ax.plot(energyList,recoEfficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  # ax.set_ylabel(r"trigger efficiency", fontsize=22)
  ax.set_ylabel(r"reconstruction efficiency", fontsize=22)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  ax.text(0.78,0.1,s=r"trig:{0}".format(triggerType),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_xscale('log')
  ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  ax.set_ylim(0,1.01)
  ax.set_xlim(14,17)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  ax.grid(True,alpha=0.5)
  l1=ax.legend(loc="upper left",fontsize=12)
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"0.98")
  # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  # ax.add_artist(l1)
  # ax.add_artist(l2)
  plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"Reco"+str(int(recoThreshold))+str(denoType)+"Efficiency.pdf",transparent=False,bbox_inches='tight')
  plt.close()

# for itrigger in triggerListSelect:
#   plotRecoEfficiency(evtList,energyBins,triggerType=itrigger,containment=True,recoThreshold=1.0,denoType="trig")
#   plotRecoEfficiency(evtList,energyBins,triggerType=itrigger,containment=False,recoThreshold=1.0,denoType="trig")
#   plotRecoEfficiency(evtList,energyBins,triggerType=itrigger,containment=True,recoThreshold=10.0,denoType="trig")
#   plotRecoEfficiency(evtList,energyBins,triggerType=itrigger,containment=False,recoThreshold=10.0,denoType="trig")
#   plotRecoEfficiency(evtList,energyBins,triggerType=itrigger,containment=True,recoThreshold=20.0,denoType="trig")
#   plotRecoEfficiency(evtList,energyBins,triggerType=itrigger,containment=False,recoThreshold=20.0,denoType="trig")
#   plotRecoEfficiency(evtList,energyBins,triggerType=itrigger,containment=True,recoThreshold=360.0,denoType="trig")
#   plotRecoEfficiency(evtList,energyBins,triggerType=itrigger,containment=False,recoThreshold=360.0,denoType="trig")


def plotRecodiff(evtList,energyBins,triggerType,containment):
  '''
  plots reco efficiency in different zenith bins
  '''
  print("plotting reco efficiency for ",triggerType)
  if containment == True:
    # evtList = containedEvents(evtList,640)
    evtList = containedEvents(evtList,410)
  colorIter = iter(colorsCustom+colorsCustom)
  # colorIter = iter(colorsList)
  sin2ZenBins = [0.0,0.5,0.822]
  energyBins = 10**np.linspace(5, 8.0, 4)
  bins = np.linspace(-60,60,242)
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
    energyList = []
    efficiencyList = []
    recoEfficiencyList = []
    ncolor = colorsCustom2[nbin]
    for ebin, ebinStart in enumerate(energyBins[:-1]):
      lowEdge_E = energyBins[ebin]
      highEdge_E = energyBins[ebin+1]
      evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
      zenithDiff = [(ievt.zenithReco-ievt.zenith)*180.0/np.pi for ievt in evtEBin if not np.isnan(ievt.zenithReco)]
      # totalEvts = len(evtEBin)
      fig = plt.figure(figsize=(8,5))
      gs = gridspec.GridSpec(nrows=1,ncols=1)
      ax = fig.add_subplot(gs[0])
      ax.hist(zenithDiff,histtype="step",bins=bins,label="tank7",lw=2.5,alpha=1)   
      ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
      ax.set_xlabel(r"diff [$^{\circ}$]", fontsize=22)
      ax.set_ylabel("count", fontsize=22)
      ax.set_xlim(-60,60)
      ax.set_yscale("log")
      ax.legend(fontsize=12)
      plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"Recodiff{:.1f}_{:.1f}.pdf".format(np.arcsin(np.sqrt(binStart))*180/np.pi,np.log10(ebinStart*10**9)),transparent=False,bbox_inches='tight')
      plt.close()

for itrigger in triggerListSelect:
  plotRecodiff(evtList,energyBins,triggerType=itrigger,containment=True)






def plotSteps(triggeredEvts,ax,legendLabel,energyScale,ncolor):
  hitBins = np.linspace(14.0,17.0,31)
  energy = [ievt.energy for ievt in triggeredEvts]
  weights = [ievt.H4aWeight for ievt in triggeredEvts]
  # weights_direct = [ievt.directWeight for ievt in triggeredEvts]
  energy = np.log10(energy)+9.0
  hist,binEdge = np.histogram(energy,hitBins,weights=[w for w in weights])
  # binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
  print("sum of simulated rate",sum(hist))
  # if "total" in legendLabel:
    # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
  # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
    # ax.set_ylim(10**-4,20)
  binCenter = (binEdge[:-1]+binEdge[1:])/2.0
  if str(energyScale) == "0.0":
    H = [h for h in hist]   
  else:
    H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
  # ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
  ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.2f} Hz".format(sum(hist)),color=ncolor,alpha=1)
  return ax,sum(hist)


def plotRecoEnergyFlux(eventList,triggerType,yscale,suffix,energyScale,containment):
  """
  plots energy flux
  """
  eventList = [ievt for ievt in eventList if (abs(getattr(ievt,triggerType)-1)<0.01 and abs(ievt.zenithReco) <= np.pi)]
  if containment == True:
    eventList = containedEvents(eventList,410)    
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(ncols=1,nrows=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  hitBins = np.linspace(14.0,17.0,31)
  totalRate = 0
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
    energyList = []
    neventList = []
    # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
    ncolor = next(colorIter)
    ax,histSum = plotSteps(evtZenBin,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),energyScale,ncolor=colorsCustom2[nbin])
    totalRate += histSum
  # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
  # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
  ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
  ax.set_yscale(yscale) 
  ax.set_ylim(10**-5,10**1)
  ax.set_xlim(14.0,17.0)
  # ax.set_xscale('log')
  ax.grid(True,alpha=0.7)
  ax.text(0.82,0.82,s=r"{0} rate:{1:.1f} Hz".format(triggerType,totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
  # ax.set_title(key,fontsize=16)
  ax.legend(fontsize=10,ncol=3,loc="lower center")
  plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+"Reco.pdf",transparent=False,bbox_inches='tight')
  plt.close()


# for itrigger in triggerList:
#   plotRecoEnergyFlux(evtList,triggerType=itrigger,yscale="log",suffix="fluxLog",energyScale=0.0,containment=False)
#   plotRecoEnergyFlux(evtList,triggerType=itrigger,yscale="log",suffix="fluxLog",energyScale=0.0,containment=True)