#!/usr/bin/env python3

import os
import glob
import subprocess
import re

import tables
import pandas as pd
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
from scipy.optimize import curve_fit

from astropy.stats import binom_conf_interval
# from statsmodels.stats.proportion import proportion_confint


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

# from weighting import GetWeight, ParticleType, PDGCode
from inclinedTriggerTools import *

# from icecube.weighting.fluxes import GaisserH4a_IT
from weighting.python.fluxes import GaisserH4a_IT
# from weighting.fluxes import GaisserH4a_IT
from weighting.python.weighting import icetop_mc_weights

# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanSeedSame/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/"
basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanFRT/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTest/"


plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

# pickleFilesPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetPickle/"
pickleFilesPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetPickle1_6/"
pickleFiles = sorted(glob.glob(pickleFilesPath+"*GenDetFiltProcUniqueCleanVEMEvts.pkl"))


evtList = []
for ipickle in pickleFiles:
  with open(ipickle,"rb") as f:
    ievtList = pickle.load(f)
    evtList += ievtList
firstEvent = evtList[0]
print("first event",dir(firstEvent))



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

# triggerList = ["HLC6_5000","tank6_5000","tank6_4000","tank6_3000","tank6_2000",
#   "tank7_5000","tank7_4000","tank7_3000","tank7_2000","tank8_5000","tank8_4000","tank8_3000","tank8_2000"]
# triggerList7 = ["HLC6_5000","tank7_5000","tank7_4000","tank7_3000","tank7_2000"]
# triggerListSelect = ["HLC6_5000","tank6_3000","tank7_3000","tank8_3000"]
# triggerListSelect = ["SMT102","SMT173","SMT183","SMT273","SMT373"]
triggerListSelect = ["SMT102","SMT273"]
# triggerListSelectDict = {"SMT102":"6HLC_5","SMT173":"7Hit_3","SMT183":"8Hit_3"}
# triggerListSelectDict = {"SMT102":"6HLC_5","SMT273":"7HGHit_3","SMT183":"8Hit_3"}
# triggerListSelectDict = {"SMT102":"ITSMT","SMT273":"IT7HG"}
# triggerListSelectDict = {"HG7_3000":"IT7HG","HLC6_5000":"ITSMT"}
# triggerListSelectDict = {"HG7_3000":"current_trig","HLC6_5000":"new_trig"}
triggerListSelectDict = {"HG7_3000":"7 tanks","HLC6_5000":"3 stations"}
# triggerListSelectDict = {"SMT102":"current trig","SMT273":"new trig"}
# triggerListFull = ["SMT102", "SMT103", "SMT104", "SMT1101", "SMT1102", "SMT1103", "SMT1104", "SMT1105", "SMT161",
#      "SMT162", "SMT163", "SMT164", "SMT165", "SMT171", "SMT172", "SMT173", "SMT174", "SMT175", "SMT181", "SMT182", "SMT183",
#       "SMT184", "SMT185", "SMT191", "SMT192", "SMT193", "SMT194", "SMT195", "SMT2101", "SMT2102", "SMT2103", "SMT2104", "SMT2105",
#        "SMT261", "SMT262", "SMT263", "SMT264", "SMT265", "SMT271", "SMT272", "SMT273","SMT274", "SMT275", "SMT281", "SMT282", "SMT283",
#         "SMT284", "SMT285", "SMT291", "SMT292", "SMT293", "SMT294", "SMT295", "SMT3101", "SMT3102", "SMT3103", "SMT3104", "SMT3105",
#          "SMT321", "SMT322", "SMT323", "SMT324", "SMT325", "SMT341", "SMT342", "SMT343", "SMT344", "SMT345", "SMT361", "SMT362",
#           "SMT363", "SMT364", "SMT365", "SMT371", "SMT372", "SMT373", "SMT374", "SMT375", "SMT381", "SMT382", "SMT383", "SMT384",
#            "SMT385","SMT391","SMT392","SMT393","SMT394","SMT395","HLC6_5000","tank6_5000","tank6_4000","tank6_3000","tank6_2000",
#   "tank7_5000","tank7_4000","tank7_3000","tank7_2000","tank8_5000","tank8_4000","tank8_3000","tank8_2000"]

# triggerListFull = ["SMT102","SMT103","SMT104","SMT163","SMT173","SMT183","SMT193","SMT263","SMT273","SMT283","SMT293","SMT323","SMT343", "SMT344", "SMT345", "SMT361", "SMT362",
#           "SMT363","SMT373","SMT383","SMT393"]


trigWindow = 10**(-6) # in ns

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
# sin2ZenBins = [0.0,0.822]

# evtList = extractEvents(hdf5NullList)
# evtListFe = extractEvents(hdf5NullListFe)
# evtListP = extractEvents(hdf5NullListP)
# print("event list before",len(evtList))
# evtList = removeFirstCore(evtList)
# print("event list after",len(evtList))
# SLCRate(evtList)
# evtList_contained = containedEvents(evtList,640)
evtList_contained = containedEvents(evtList,410)
total_events = sum([1 for ievt in evtList])
print("total events",total_events)

# weightsI3 = [ievt.H4aWeight for ievt in evtList]
# weightsPy1 = [ievt.H4aWeight2 for ievt in evtList]
# weightsPy2 = weightCalc(hdf5NullList)
# weightsPy3 = weightCalc1(hdf5NullList,nfilesP,nfilesHe,nfilesO,nfilesFe)
# triggered_events = sum([ievt.ITSMTTriggered for ievt in evtList])
# print("trigger rate",total_events,triggered_events,triggered_events/total_events)
def plotWeightHist(weightsi3,weightspy2,weightspy3):
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  print("length",len(weightsi3),min(weightsi3),max(weightsi3))
  print("length",len(weightspy2),min(weightspy2),max(weightspy2))
  print("length",len(weightspy3),min(weightspy3),max(weightspy3))
  print("ratio",len(weightspy3),min(weightspy3)/min(weightsi3),max(weightspy3)/max(weightsi3))
  # ax.hist(weightsi3,histtype="step")
  ax.hist(weightspy2,histtype="step")
  # ax.hist(weightspy3,histtype="step")
  plt.savefig(plotFolder+"/weightCompare3.pdf",transparent=False,bbox_inches='tight')
  plt.close()

# plotWeightHist(weightsI3,weightsPy2,weightsPy3)


# # weights = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zenith,energy,ptype)
# # weithtsOne = GetWeight().getWeight(nfilesP,nfilesHe,nfilesO,nfilesFe,zenith[-1],energy[-1],ptype[-1])
# # print("weightCompare",weights[-1],weithtsOne,len(weights))
# # adjustedWeights = []



energyBins = 10**np.linspace(5, 8.0, 31)
energyBinsShort = 10**np.linspace(5, 8.0, 7)
# energyBins = 10**np.linspace(5, 8.0, 7)
energyBinslgE = np.linspace(5.0,8.9,4000)
energyBinCenter = [5.1,6.1,7.1,8.1]
print("energy bins",energyBins)


def sigmoid(x, a,b,c):
  """returns a sigmoid function """
  # return [1 / (1 + np.exp(-b*(ix-a)))+c for ix in x]
  x = np.asarray(x)
  return 1 / (1 + np.exp(-b*(x-a)))+c



# def sigmoid(x, a,b):
#   """returns a sigmoid function """
#   # return [1 / (1 + np.exp(-b*(ix-a)))+c for ix in x]
#   x = np.asarray(x)
#   # print("debug",x)
#   print("debug a",a)
#   print("debug b",b)
#   return 1 / (1 + np.exp(-b*(x-a)))

def poissonError(n):
  return np.sqrt(n)

def binomial_proportion(nsel, ntot, coverage=0.68):
    """
    Calculate a binomial proportion (e.g. efficiency of a selection) and its confidence interval.

    Parameters
    ----------
    nsel: array-like
      Number of selected events.
    ntot: array-like
      Total number of events.
    coverage: float (optional)
      Requested fractional coverage of interval (default: 0.68).

    Returns
    -------
    p: array of dtype float
      Binomial fraction.
    dpl: array of dtype float
      Lower uncertainty delta (p - pLow).
    dpu: array of dtype float
      Upper uncertainty delta (pUp - p).
    """

    from scipy.stats import norm

    z = norm().ppf(0.5 + 0.5 * coverage)
    z2 = z * z
    p = np.asarray(nsel, dtype=float) / ntot
    div = 1.0 + z2 / ntot
    pm = (p + z2 / (2 * ntot))
    dp = z * np.sqrt(p * (1.0 - p) / ntot + z2 / (4 * ntot * ntot))
    pl = (pm - dp) / div
    pu = (pm + dp) / div

    return p, p - pl, pu - p

####################################################
# # Short tools to help define Wilson confidence intervals for binomial trials
# # The z-parameter gives the equivalent z-score for a Gaussian distribution
# # z = 1.00, 68%
# # z = 1.44, 85%
# # z = 1.96, 95%


def WilsonMean(nPass, nFail, z=1.0):
    if nPass < 0 or nFail < 0:
        return 0

    return (nPass + z ** 2 / 2.0) / float(nPass + nFail + z ** 2)


def WilsonError(nPass, nFail, z=1.0):
    if nPass < 0 or nFail < 0:
        return 0

    if nPass + nFail == 0:
        return 0.5

    term = z / float(nPass + nFail + z ** 2)

    return term * (nPass * nFail / float(nPass + nFail) + z ** 2 / 4.0) ** 0.5

###################################################



from iminuit import Minuit

class MinuitFit(object):
  """docstring for MinuitFit"""
  def __init__(self,x,y,yerr,fitFunc):
    super(MinuitFit, self).__init__()
    self.x = x
    self.y = y
    self.yerr = yerr
    self.fitFunc = fitFunc

  def AccuracyMetric(self,expected,observed,error,n_parameter):
    # return np.sum(((expected - observed)/error)**2) / (len(observed)-n_parameter)
    sum_metric = 0
    for e,o,er in zip(expected,observed,error):
      if er > 0:
        sum_metric += ((e - o)/er)**2 / len(observed)
      else:
        sum_metric += ((e - o)/0.00001)**2 / len(observed)
    return sum_metric

  def LossFcn(self,a,b,c):
    # expected = self.fitFunc(self.x,*par)
    expected = np.asarray(self.fitFunc(self.x,a,b,c))
    observed = np.asarray(self.y)
    error = []
    if len(np.asarray(self.yerr).shape) == 1:
      error = self.yerr
    elif len(np.asarray(self.yerr).shape) == 2:
      for ielt in range(len(self.x)):
        if expected[ielt] >= observed[ielt]:
          error.append(self.yerr[1][ielt])
        else:
          error.append(self.yerr[0][ielt])
    n_parameter = 3
    error = np.asarray(error)
    return self.AccuracyMetric(expected,observed,error,n_parameter)

  def minimization(self):
    # m = Minuit(self.LossFcn,a=14.0,error_a = 0.2,b=0.5,error_b = 0.2,limit_a=(1.0,18.0),limit_b=(0.01,0.99)
    #   ,c=0.1,error_c=0.2,limit_c=(-100.0,100.0),errordef=1,pedantic=False)
    m = Minuit(self.LossFcn,a=14.0,error_a = 0.2,b=7.0,error_b = 0.2,c=0.1,error_c=0.2,
      limit_a=(14.0,18.0),limit_b=(4.0,12.0),errordef=1,pedantic=False)
    # m.limits['a'] = (13.0,14.0)
    # m.limits['b'] = (0.1,0.9)
    # ,c=0.0,error_c = 0.2
    # m = Minuit.from_array_func(self.LossFcn,init = [14.0,0.5,0.1],limit=[[13.0,15.0],[0.1,0.9]], name=["a","b"], error=[0.01,0.01], errordef=1.,pedantic=False)
    # m = Minuit.from_array_func(self.LossFcn,init = [14.0,0.5,0.1],name=["a","b","c"],errordef=1.,pedantic=False)
    # m = Minuit(self.LossFcn,name=name,errordef=1.,pedantic=False)
    m.migrad() #find minimum of loss function
    # loss = self.fitFunc(m.np_values())
    loss = None
    print("parameters",m.parameters,m.np_values())
    print("covariance", m.np_covariance())
    return m, loss





def plotTrigEfficiencyPure(eventList,energyBins,triggerType,containment,sigFit,errorType):
  '''
  plots trigger efficiency in different zenith bins
  '''
  print("plotting trigger efficiency for ",triggerType)
  eventListP = [ievt for ievt in eventList if str(ievt.CRType) == str(2212)]
  eventListFe = [ievt for ievt in eventList if str(ievt.CRType) == str(1000260560)]
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # colorIter = iter(colorsList)
  # sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
  sin2ZenBins = [0.0,0.5,0.822]
  # sin2ZenBins = sin2ZenBins[-4:]
  colorListp = ["maroon","firebrick","indianred","lightcoral"]
  colorListFe = ["darkslategray","teal","darkturquoise","skyblue"]
  pMap = {"p":eventListP,"Fe":eventListFe}
  cMap = {"p":colorListp,"Fe":colorListFe}
  for iprimary in pMap.keys():
    eventList = pMap[iprimary]
    if containment == True:
      # evtList = containedEvents(evtList,640)
      # eventList = containedEvents(eventList,410)
      eventList = containedEvents(eventList,800)
    for nbin, binStart in enumerate(sin2ZenBins[:-1]):
      lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
      highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
      evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
      energyList = []
      efficiencyList = []
      errorListLow = []
      errorListHigh = []
      ncolor = cMap[iprimary][nbin]
      for ebin, ebinStart in enumerate(energyBins[:-1]):
        lowEdge_E = energyBins[ebin]
        highEdge_E = energyBins[ebin+1]
        evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
        # totalEvts = len(evtEBin)
        weights = [ievt.H4aWeight for ievt in evtEBin]
        totalEvts = len(evtEBin)
        # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
        triggerList = [ievt for ievt in evtEBin if abs(getattr(ievt,triggerType)-1)<0.01]
        trigEff = triggerEfficiency(len(triggerList),totalEvts)
        efficiencyList.append(trigEff)
        if errorType == "poisson":            
          if totalEvts != 0:
            poissonErr = np.sqrt(len(triggerList))/totalEvts
          else:
            poissonErr = 0.00001
          errorListLow.append(poissonErr)
          errorListHigh.append(poissonErr)
        elif errorType == "wilson":
          # wilsonErr = binom_conf_interval(k=len(triggerList), n=totalEvts, confidence_level=0.68269, interval='wilson')
          if totalEvts < 1:
            errorListLow.append(0.000001)
            errorListHigh.append(0.000001)
          else:
            # wilsonErr = binom_conf_interval(k=len(triggerList), n=totalEvts, interval='wilson')
            wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68)
            errorListLow.append(wilsonErr[1])
            errorListHigh.append(wilsonErr[2])
        energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))
      if sigFit == True:
        ###########################################scipy fit#################
        # popt, pcov = curve_fit(sigmoid, energyList[:-2], efficiencyList[:-2], p0=[14, 0.5,0])
        # xfit = np.linspace(14,17,1000)
        # yfit = sigmoid(xfit,*popt)
        # ax.plot(xfit,yfit,ls="-",c=ncolor,alpha=1.0)
        # # ax.plot(energyList,efficiencyList,".",lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$ {2:}".format(
        # #   np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,iprimary),alpha=1)
        # ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$ {2:}".format(
        #   np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,iprimary),alpha=1)
        # ###########################################minuit fit################
        xfit = np.linspace(14,17,1000)
        yerr=np.asarray([errorListLow[:-8],errorListHigh[:-8]])
        fit = MinuitFit(x=energyList[:-8],y=efficiencyList[:-8],yerr=yerr,fitFunc=sigmoid)
        m,loss = fit.minimization()
        print("m values",*m.np_values())
        ax.plot(xfit,sigmoid(xfit,*m.np_values()),ls="-",lw=2.5,c=ncolor,alpha=1.0)
        ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$ {2:}".format(
          np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,iprimary),alpha=1)
      else:
        ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$ {2:}".format(
          np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,iprimary),alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  ax.set_ylabel(r"trigger efficiency", fontsize=22)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  ax.text(0.78,0.1,s=r"trig:{0}".format(triggerListSelectDict[triggerType]),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_xscale('log')
  ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  ax.set_ylim(0,1.05)
  ax.set_xlim(14,17)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  ax.grid(True,alpha=0.5)
  l1=ax.legend(loc='center right',fontsize=10)
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"0.98")
  # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  # ax.add_artist(l1)
  # ax.add_artist(l2)
  plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"PureEfficiency.pdf",transparent=False,bbox_inches='tight')
  plt.close()

# for itrigger in triggerListSelectDict.keys():
#   plotTrigEfficiencyPure(evtList,energyBins,triggerType=itrigger,containment=True,sigFit=False,errorType="poisson")
# plotTrigEfficiencyPure(evtList,energyBins,triggerType="tank7_3000",containment=True,sigFit=True,errorType="wilson")


def nHitsPerEvent(eventList,bins,hitType,triggerType,containment):
  if containment == True:
    # evtList = containedEvents(evtList,640)
    # eventList = containedEvents(eventList,410)
    eventList = containedEvents(eventList,800)
  nHits = [getattr(ievt,hitType) for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
  meanHit = np.mean(nHits)
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  bins = np.linspace(-1,max(nHits),max(nHits)+2)
  ax.hist(nHits,bins=bins,histtype="step",label="total",lw=1.5)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel("n_hits per event", fontsize=22)
  ax.set_ylabel("count", fontsize=22)
  ax.set_yscale("log")
  ax.grid(True,alpha=0.4)
  pulseDict = {"nSLC":"SLC tank","nHLC":"HLC tank","nHLCVEM":"HLC VEM","nSLCVEM":"SLC VEM"}
  ax.text(0.82,0.88,s="mean {0} hits\n{1:.1f} per event".format(pulseDict[hitType],meanHit),size=15,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.legend(loc="upper left",fontsize=12)
  # ax.legend(fontsize=12)
  # ax.hist(weightspy3,histtype="step")
  plt.savefig(plotFolder+"/hits"+str(hitType)+str(triggerType)+str(containment)+".pdf",transparent=False,bbox_inches='tight')
  plt.close()

# nHitsPerEvent(evtList,bins=np.linspace(-1,200,202),hitType="nHLCVEM",triggerType="tank7_3000",containment=True)
# nHitsPerEvent(evtList,bins=np.linspace(-1,60,62),hitType="nSLCVEM",triggerType="tank7_3000",containment=True)
# nHitsPerEvent(evtList,bins=np.linspace(-1,200,202),hitType="nHLC",triggerType="tank7_3000",containment=True)
# nHitsPerEvent(evtList,bins=np.linspace(-1,60,62),hitType="nSLC",triggerType="tank7_3000",containment=True)
# nHitsPerEvent(evtList,bins=np.linspace(-1,200,202),hitType="nHLCVEM",triggerType="tank7_3000",containment=False)
# nHitsPerEvent(evtList,bins=np.linspace(-1,60,62),hitType="nSLCVEM",triggerType="tank7_3000",containment=False)
# nHitsPerEvent(evtList,bins=np.linspace(-1,200,202),hitType="nHLC",triggerType="tank7_3000",containment=False)
# nHitsPerEvent(evtList,bins=np.linspace(-1,60,62),hitType="nSLC",triggerType="tank7_3000",containment=False)

# nHitsPerEvent(evtList,bins=np.linspace(-1,200,202),hitType="nHLCVEM",triggerType="HLC6_5000",containment=True)
# nHitsPerEvent(evtList,bins=np.linspace(-1,60,62),hitType="nSLCVEM",triggerType="HLC6_5000",containment=True)
# nHitsPerEvent(evtList,bins=np.linspace(-1,200,202),hitType="nHLC",triggerType="HLC6_5000",containment=True)
# nHitsPerEvent(evtList,bins=np.linspace(-1,60,62),hitType="nSLC",triggerType="HLC6_5000",containment=True)
# nHitsPerEvent(evtList,bins=np.linspace(-1,200,202),hitType="nHLCVEM",triggerType="HLC6_5000",containment=False)
# nHitsPerEvent(evtList,bins=np.linspace(-1,60,62),hitType="nSLCVEM",triggerType="HLC6_5000",containment=False)
# nHitsPerEvent(evtList,bins=np.linspace(-1,200,202),hitType="nHLC",triggerType="HLC6_5000",containment=False)
# nHitsPerEvent(evtList,bins=np.linspace(-1,60,62),hitType="nSLC",triggerType="HLC6_5000",containment=False)



# plotRadiusEnergy(energyBinslgE)

# def plotCoreScatter_(x,y,suffix,title):
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   # ax.scatter(x,y,s=10,alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
#   ax.set_xlabel(r"x [m]", fontsize=24)
#   ax.set_ylabel(r"y [m]", fontsize=24)
#   xCirc,yCirc = getCircle(800)
#   ax.plot(xCirc,yCirc,'-',c="purple",lw=3.0,label="r = 800 m")
#   xCirc,yCirc = getCircle(1100)
#   ax.plot(xCirc,yCirc,'-',c="blue",lw=3.0,label="r = 1100 m")
#   xCirc,yCirc = getCircle(1700)
#   ax.plot(xCirc,yCirc,'-',c="orange",lw=3.0,label="r = 1700 m")
#   # xCirc,yCirc = getCircle(2600)
#   # ax.plot(xCirc,yCirc,'-',c="yellow",lw=3.0,label="r = 2600 m")
#   ax.scatter(x,y,s=10,alpha=1)
#   # ax.set_xlim(0,100)
#   # ax.set_ylim(0,100)
#   ax.grid(True,alpha=0.2)
#   ax.set_title(title,fontsize=16)
#   ax.set_aspect("equal")
#   plt.legend(fontsize=12)
#   plt.savefig(plotFolder+"/coreScatter"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # def plotCoreScatter(hdfFileList):
# #   x,y = getCore(hdfFileList)
# #   zenList,ptypeList,energyList = getZenithTypeEnergy(hdfFileList)
# #   plotCoreScatter_(x,y,"all_shower","all energy")
# #   energyBins = 10**(np.linspace(5.0,8.0,31))
# #   print("energyBins",energyBins)
# #   for n,nEnergy in enumerate(energyBins[:-1]):
# #     xInBin = []
# #     yInBin = []
# #     for ix,iy,ienergy in zip(x,y,energyList):
# #       if ienergy >= energyBins[n] and ienergy < energyBins[n+1]:
# #         xInBin.append(ix) 
# #         yInBin.append(iy)
# #     print("nBins",n,energyBins[n],energyBins[n+1])
# #     plotCoreScatter_(xInBin,yInBin,r"{:.1f}".format(np.log10(nEnergy)),r"lg(E[GeV]):{0:.1f}-{1:.1f}".format(np.log10(energyBins[n]),np.log10(energyBins[n+1])))

# def plotCoreScatterEnergy(evtList,energyLow,energyHigh,filtKey):
#   x = [ievt.coreX for ievt in evtList]
#   y = [ievt.coreY for ievt in evtList]
#   zenithBins = [np.arcsin(np.sqrt(i)) for i in np.linspace(0.0,1.0,11)][:-3]
#   zenithBins.append(np.deg2rad(65))
#   print("zenith bins",[np.sin(i)**2 for i in zenithBins])
#   for n,nZenith in enumerate(zenithBins[:-1]):
#     xInBin = []
#     yInBin = []
#     for ix,iy,izen,ienergy in zip(x,y,zenList,energyList):
#       if ienergy >= 10**energyLow and ienergy < 10**energyHigh and izen >= zenithBins[n] and izen < zenithBins[n+1]:
#         xInBin.append(ix) 
#         yInBin.append(iy)
#     plotCoreScatter_(xInBin,yInBin,r"energy_{:.1f}Zen{:.1f}_filt{}".format(energyLow,np.sin(nZenith)**2,str(filtKey)),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(np.rad2deg(zenithBins[n]),np.rad2deg(zenithBins[n+1])))

# # plotCoreScatterEnergy(evtList,6.0,6.1,"None")
# # plotCoreScatterEnergy([ievt for ievt in evtList if abs(ievt.ITSMTTriggered-1)<0.01] ,6.0,6.1,"sta3")
# # plotCoreScatterEnergy(evtList,7.0,7.1,"None")
# # plotCoreScatterEnergy([ievt for ievt in evtList if abs(ievt.ITSMTTriggered-1)<0.01],7.0,7.1,"sta3")
# # plotCoreScatterEnergy(evtList,5.9,6.0,"None")
# # plotCoreScatterEnergy([ievt for ievt in evtList if abs(ievt.ITSMTTriggered-1)<0.01],6.9,7.0,"sta3")


# def plotCoreScatter_(x,y,suffix,title):
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   # ax.scatter(x,y,s=10,alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
#   ax.set_xlabel(r"x [m]", fontsize=24)
#   ax.set_ylabel(r"y [m]", fontsize=24)
#   xCirc,yCirc = getCircle(800)
#   ax.plot(xCirc,yCirc,'-',c="purple",lw=3.0,label="r = 800 m")
#   xCirc,yCirc = getCircle(1100)
#   ax.plot(xCirc,yCirc,'-',c="blue",lw=3.0,label="r = 1100 m")
#   xCirc,yCirc = getCircle(1700)
#   ax.plot(xCirc,yCirc,'-',c="orange",lw=3.0,label="r = 1700 m")
#   # xCirc,yCirc = getCircle(2600)
#   # ax.plot(xCirc,yCirc,'-',c="yellow",lw=3.0,label="r = 2600 m")
#   ax.scatter(x,y,s=0.05,alpha=1)
#   # ax.set_xlim(0,100)
#   # ax.set_ylim(0,100)
#   ax.grid(True,alpha=0.2)
#   ax.set_title(title,fontsize=30)
#   ax.set_aspect("equal")
#   plt.legend(fontsize=12)
#   plt.savefig(plotFolder+"/coreScatter"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()


# def plotScatterCore(evtList,triggerType):
#   '''
#   plots core distance of triggered shower
#   '''
#   energyBins = 10**(np.linspace(5.0,8.0,4))
#   evtList = selectTriggered(evtList,triggerType)
#   for ebin, ebinStart in enumerate(energyBins[:-1]):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     colorIter = iter(colorsList)
#     lowEdge_E = energyBins[ebin]
#     highEdge_E = energyBins[ebin+1]
#     evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
#     ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((energyBins[ebin])*10**9),np.log10((energyBins[ebin+1])*10**9)))
#     for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#       # colorIter = iter(colorsCustom)
#       lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#       highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#       evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
#       distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
#       x = [ievt.coreX for ievt in evtZenBin]
#       y = [ievt.coreY for ievt in evtZenBin]
#       if len(distanceList)>2:
#         xbins = np.linspace(min(distanceList),max(distanceList),200)
#       else:
#         xbins = np.linspace(0,1000,200)
#       ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
#         r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),lw=2.5,alpha=1)
#       plotCoreScatter_(x,y,r"energy_{:.1f}Zen{:.1f}Trig{}".format(np.log10(lowEdge_E),np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180/np.pi,triggerType),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180/np.pi))
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#     ax.set_xlabel(r"core distance [m]", fontsize=22)
#     ax.set_ylabel(r"count", fontsize=22)
#     # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#     # ax.set_xscale('log')
#     ax.set_ylim(0,1300)
#     # ax.set_ylim(0,600)
#     ax.set_xlim(0,1710)
#     # ax.yaxis.set_minor_locator(MultipleLocator(100))
#     # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#     ax.grid(True,alpha=0.6)
#     ax.legend(fontsize=12,loc='upper left')
#     plt.savefig(plotFolder+"/distanceScattTrig"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(energyBins[ebin])+9),transparent=False,bbox_inches='tight')
#     plt.close()
# plotScatterCore(evtList,"None")
# plotScatterCore(evtList,"sta3")
# plotScatterCore(evtList,"sta1")

# def plot2dHistCoreGiven(evtList,triggerType,energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-800,800],ylim=[-800,800]):
#   '''
#   plots core distance of triggered shower as a 2d histogram
#   energy Lim given in eV.
#   '''
#   energyBins = 10**(np.linspace(5.0,8.0,4))
#   evtList = selectTriggered(evtList,triggerType)
#   fig = plt.figure(figsize=(8,8))
#   # fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsList)
#   lowEdge_E = energyLim[0]*10**(-9)
#   highEdge_E = energyLim[1]*10**(-9)
#   evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
#   ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}, trig:{2}".format(np.log10((lowEdge_E)*10**9),np.log10((highEdge_E)*10**9),triggerType))
#   # colorIter = iter(colorsCustom)
#   lowEdge = zenLim[0]
#   highEdge = zenLim[1]
#   evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
#   distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
#   # idBug = [(ievt.runID,ievt.eventID,ievt.coreX,ievt.coreY) for ievt in evtZenBin if abs(ievt.coreX-444)<0.4 and abs(ievt.coreY-0)<1]
#   # print("xyBug",idBug)
#   x = [ievt.coreX for ievt in evtZenBin]
#   y = [ievt.coreY for ievt in evtZenBin]
#   xbins = np.linspace(min(distanceList),max(distanceList),200)
#   # ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
#   #   r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(zenLim[0]*180.0/np.pi,zenLim[1]*180.0/np.pi),lw=2.5,alpha=0.4)
#   # plotCoreScatter_(x,y,r"energy_{:.1f}Zen{:.1f}".format(lowEdge_E,zenLim[0]*180.0/np.pi),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(zenLim[0]*180.0/np.pi,zenLim[0]*180.0/np.pi))
#   counts, xedges, yedges, im = ax.hist2d(x,y,bins=100,norm=mpl.colors.LogNorm())
#   cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
#   cbar.set_label('# of cores',fontsize=18)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"core distance x[m]", fontsize=22)
#   ax.set_ylabel(r"core distance y[m]", fontsize=22)
#   ax.set_aspect("equal")
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   # ax.set_xscale('log')
#   # ax.set_ylim(-810,810)
#   # ax.set_xlim(-810,810)
#   ax.set_ylim(ylim[0],ylim[1])
#   ax.set_xlim(xlim[0],xlim[1])
#   # ax.yaxis.set_minor_locator(MultipleLocator(100))
#   # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#   ax.grid(True,alpha=0.6)
#   # ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/distance2dHist"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(lowEdge_E)+9),transparent=False,bbox_inches='tight')
#   plt.close()


# # plot2dHistCoreGiven(evtList,"none",energyLim=[10**14.0,10**15.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-800,800],ylim=[-800,800])
# # plot2dHistCoreGiven(evtList,"none",energyLim=[10**15.0,10**16.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-1100,1100],ylim=[-1100,1100])
# # plot2dHistCoreGiven(evtList,"none",energyLim=[10**16.0,10**17.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-1700,1700],ylim=[-1700,1700])

# # plot2dHistCoreGiven(evtList,"sta3",energyLim=[10**14.0,10**15.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-800,800],ylim=[-800,800])
# # plot2dHistCoreGiven(evtList,"sta3",energyLim=[10**15.0,10**16.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-1100,1100],ylim=[-1100,1100])
# # plot2dHistCoreGiven(evtList,"sta3",energyLim=[10**16.0,10**17.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-1700,1700],ylim=[-1700,1700])

# # plot2dHistCoreGiven(evtList,"sta1",energyLim=[10**14.0,10**15.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-800,800],ylim=[-800,800])
# # plot2dHistCoreGiven(evtList,"sta1",energyLim=[10**15.0,10**16.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-1100,1100],ylim=[-1100,1100])
# # plot2dHistCoreGiven(evtList,"sta1",energyLim=[10**16.0,10**17.0],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],xlim=[-1700,1700],ylim=[-1700,1700])


# def plotScatterCoreGiven(evtList,triggerType,energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[400,500]):
#   '''
#   plots core distance of triggered shower
#   '''
#   energyBins = 10**(np.linspace(5.0,8.0,7))
#   # evtList = selectTriggered(evtList,triggerType)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsList)
#   lowEdge_E = energyLim[0]
#   highEdge_E = energyLim[1]
#   evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
#   ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((lowEdge_E)*10**9),np.log10((highEdge_E)*10**9)))
#   # colorIter = iter(colorsCustom)
#   lowEdge = zenLim[0]
#   highEdge = zenLim[1]
#   evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
#   distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
#   x = [ievt.coreX for ievt in evtZenBin]
#   y = [ievt.coreY for ievt in evtZenBin]
#   xbins = np.linspace(min(distanceList),max(distanceList),200)
#   ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
#     r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(zenLim[0]*180.0/np.pi,zenLim[1]*180.0/np.pi),lw=2.5,alpha=0.4)
#   plotCoreScatter_(x,y,r"energy_{:.1f}Zen{:.1f}".format(np.log10(lowEdge_E),zenLim[0]*180.0/np.pi),r"$\theta^{{\circ}}$:{0:.1f}-{1:.1f}".format(zenLim[0]*180.0/np.pi,zenLim[0]*180.0/np.pi))
#   # ax.hist2d(x,y,bins=100)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"core distance [m]", fontsize=22)
#   ax.set_ylabel(r"count", fontsize=22)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   # ax.set_xscale('log')
#   # ax.set_ylim(0,1.05)
#   # ax.set_xlim(14,17)
#   # ax.yaxis.set_minor_locator(MultipleLocator(100))
#   # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#   ax.grid(True,alpha=0.6)
#   ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/distanceScattTrig"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(lowEdge_E)+9),transparent=False,bbox_inches='tight')
#   plt.close()
# # plotScatterCoreGiven(evtList,"sta3",energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[400,500])







# def plotDeltaT(eventList,yscale,suffix,energyScale,triggerType):
#   """
#   plots energy flux
#   """
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   hitBins = np.linspace(14.0,17.0,31)
#   totalRate = 0
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#     highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#     evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
#     energyList = []
#     neventList = []
#     # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#     ax,histSum = plotSteps(evtZenBin,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),energyScale)
#     totalRate += histSum
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20)
#   ax.set_yscale(yscale) 
#   ax.set_ylim(10**-5,10**1)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.4)
#   ax.text(0.82,0.82,s=r"{0} rate:{1:.1f} Hz".format(triggerType,totalRate),size=15,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#   # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()


def plotEnergyFlux(eventList,triggerType,yscale,suffix,energyScale,containment):
  """
  plots energy flux
  """
  eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-0)>0.01]
  if containment == True:
    # eventList = containedEvents(eventList,410)
    # eventList = containedEvents(eventList,410)
    eventList = containedEvents(eventList,800)
  # eventList = [ievt for ievt in eventList if 10**6.9 <= ievt.energy < 10**8.0]    
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(ncols=1,nrows=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # hitBins = np.linspace(14.0,17.0,31)
  hitBins = np.linspace(14.0,16.8,29)
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
  ax.set_ylim(10**-5,10**0.0)
  # ax.set_xlim(14.0,17.0)
  ax.set_xlim(14.0,16.8)
  # ax.set_xscale('log')
  ax.grid(True,alpha=0.7)
  ax.text(0.82,0.82,s=r"{0} rate:{1:.1f} Hz".format(triggerListSelectDict[triggerType],totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
  # ax.set_title(key,fontsize=16)
  ax.legend(fontsize=10,ncol=3,loc="lower center")
  plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+".pdf",transparent=False,bbox_inches='tight')
  plt.close()

# def plotEnergyFluxComponentTest(eventList,triggerType,yscale,suffix,energyScale,containment):
#   """
#   plots energy flux
#   """
#   eventList = [ievt for ievt in eventList if abs(getattr(ievt,triggerType)-1)<0.01]
#   if containment == True:
#     eventList = containedEvents(eventList,410)    
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   hitBins = np.linspace(14.0,17.0,31)
#   totalRate = 0
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#     highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#     evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
#     energyList = []
#     neventList = []
#     # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#     ncolor = next(colorIter)
#     ax,histSum = plotSteps(evtZenBin,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),energyScale,ncolor=colorsCustom2[nbin])
#     totalRate += histSum
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
#   ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
#   ax.set_yscale(yscale) 
#   ax.set_ylim(10**-5,10**1)
#   ax.set_xlim(14.0,17.0)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.7)
#   ax.text(0.82,0.82,s=r"{0} rate:{1:.1f} Hz".format(triggerType,totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#   # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+"ComponentTest.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# def plotPureEnergyFlux(eventListP,eventListFe,triggerType,yscale,suffix,energyScale,containment):
#   """
#   plots energy flux
#   """   
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   hitBins = np.linspace(14.0,17.0,31)
#   totalRate = 0
#   sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
#   sin2ZenBins = sin2ZenBins[-4:]
#   colorListp = ["maroon","brown","indianred","lightcoral"]
#   colorListFe = ["darkslategray","teal","cyan","skyblue"]
#   pMap = {"p":eventListP,"Fe":eventListFe}
#   cMap = {"p":colorListp,"Fe":colorListFe}
#   # for eventList in [eventListP,eventListFe]:
#   for iprimary in pMap.keys():
#     eventList = [ievt for ievt in pMap[iprimary] if abs(getattr(ievt,triggerType)-1)<0.01]
#     if containment == True:
#       eventList = containedEvents(eventList,410)    
#     for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#       lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#       highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#       evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
#       energyList = []
#       neventList = []
#       # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#       ncolor = next(colorIter)
#       ax,histSum = plotSteps(evtZenBin,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$ {2:}".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,
#         np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,iprimary),energyScale,ncolor=cMap[iprimary][nbin])
#       totalRate += histSum
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
#   ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
#   ax.set_yscale(yscale) 
#   ax.set_ylim(10**-5,10**1)
#   ax.set_xlim(14.0,17.0)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.7)
#   ax.text(0.82,0.82,s=r"{0} rate:{1:.1f} Hz".format(triggerType,totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#   # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"cont"+str(containment)+"Pure.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

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
  # ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.6f} Hz".format(sum(hist)),color=ncolor,alpha=1)
  ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.3f} Hz".format(sum(hist)),color=ncolor,alpha=1)
  return ax,sum(hist)


for itrigger in triggerListSelectDict.keys():
  # plotEnergyFlux(evtList,triggerType=itrigger,yscale="log",suffix="fluxLog",energyScale=0.0,containment=False)
  plotEnergyFlux(evtList,triggerType=itrigger,yscale="log",suffix="fluxLog",energyScale=0.0,containment=True)
# for itrigger in triggerListFull:
#   plotPureEnergyFlux(evtListP,evtListFe,triggerType=itrigger,yscale="log",suffix="fluxLog",energyScale=0.0,containment=False)
#   plotPureEnergyFlux(evtListP,evtListFe,triggerType=itrigger,yscale="log",suffix="fluxLog",energyScale=0.0,containment=True)


# def plotInclinedEnergyFlux(eventList,triggerList,yscale,suffix,energyScale,containment):
#   """
#   plots energy flux
#   """
#   if containment == True:
#     eventList = containedEvents(eventList,410)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   hitBins = np.linspace(14.0,17.0,31)
#   totalRate = 0
#   lowEdge = np.arcsin(np.sqrt(sin2ZenBins[-2]))
#   highEdge = np.arcsin(np.sqrt(sin2ZenBins[-1]))
#   for ntrig,itrigger in enumerate(triggerList):
#     evtZenBin = [ievt for ievt in eventList if (lowEdge <= ievt.zenith < highEdge and abs(getattr(ievt,itrigger)-1)<0.01)]
#     energyList = []
#     neventList = []
#     # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#     ncolor = multiple_color[ntrig]
#     ax,histSum = plotSteps(evtZenBin,ax,str(itrigger),energyScale,ncolor=ncolor)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
#   ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
#   ax.set_yscale(yscale) 
#   ax.set_ylim(10**-5,10**(-1))
#   ax.set_xlim(14,17)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.7)
#   ax.text(0.78,0.1,s=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[-2]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[-1]))*180.0/np.pi),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#   # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/energySpecTrigInclined"+str(suffix)+"scale"+str(energyScale)+"Cont"+str(containment)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# plotInclinedEnergyFlux(evtList,triggerList,"log","fluxLog",0.0,containment=False)
# plotInclinedEnergyFlux(evtList,triggerList,"log","fluxLog",0.0,containment=True)

# def plotInclinedEnergyFlux7(eventList,triggerList,yscale,suffix,energyScale,containment=True):
#   """
#   plots energy flux
#   """
#   if containment == True:
#     eventList = containedEvents(eventList,410)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   hitBins = np.linspace(14.0,17.0,31)
#   totalRate = 0
#   lowEdge = np.arcsin(np.sqrt(sin2ZenBins[-2]))
#   highEdge = np.arcsin(np.sqrt(sin2ZenBins[-1]))
#   for ntrig,itrigger in enumerate(triggerList):
#     evtZenBin = [ievt for ievt in eventList if (lowEdge <= ievt.zenith < highEdge and abs(getattr(ievt,itrigger)-1)<0.01)]
#     energyList = []
#     neventList = []
#     # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#     ncolor = multiple_color[ntrig]
#     ax,histSum = plotSteps(evtZenBin,ax,str(itrigger),energyScale,ncolor=ncolor)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
#   ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
#   ax.set_yscale(yscale) 
#   ax.set_ylim(10**-5,10**(-1))
#   ax.set_xlim(14,17)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.7)
#   ax.text(0.78,0.1,s=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[-2]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[-1]))*180.0/np.pi),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#   # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/energySpecTrigInclined"+str(suffix)+"scale"+str(energyScale)+"7"+"Cont"+str(containment)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# plotInclinedEnergyFlux7(evtList,triggerList7,"log","fluxLog",0.0,containment=False)
# plotInclinedEnergyFlux7(evtList,triggerList7,"log","fluxLog",0.0,containment=True)

def plotInclinedEnergyFluxSelect(eventList,triggerListDict,yscale,suffix,energyScale,containment):
  """
  plots energy flux
  """
  triggerList = triggerListDict.keys()
  if containment == True:
    # eventList = containedEvents(eventList,410)
    eventList = containedEvents(eventList,800)
  # eventList = [ievt for ievt in eventList if 10**6.9 <= ievt.energy < 10**8.0]
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(ncols=1,nrows=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # hitBins = np.linspace(14.0,17.0,31)
  hitBins = np.linspace(14.0,16.8,29)
  totalRate = 0
  lowEdge = np.arcsin(np.sqrt(sin2ZenBins[-2]))
  highEdge = np.arcsin(np.sqrt(sin2ZenBins[-1]))
  for ntrig,itrigger in enumerate(triggerList):
    evtZenBin = [ievt for ievt in eventList if (lowEdge <= ievt.zenith < highEdge and abs(getattr(ievt,itrigger)-1)<0.01)]
    energyList = []
    neventList = []
    # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
    ncolor = multiple_color[ntrig]
    ax,histSum = plotSteps(evtZenBin,ax,str(triggerListDict[itrigger]),energyScale,ncolor=ncolor)
  # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
  # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
  ax.set_ylabel(r"rate [Hz]".format(energyScale), fontsize=20)
  ax.set_yscale(yscale) 
  ax.set_ylim(10**-4.6,10**(-1.8))
  # ax.set_ylim(10**-5,10**(-2))
  # ax.set_ylim(6*10**-5,4*10**(-3))
  # ax.set_xlim(14,17)
  ax.set_xlim(14.0,16.8)
  # ax.set_xlim(15.8,17)
  # ax.set_xscale('log')
  ax.grid(True,alpha=0.7)
  ax.text(0.78,0.1,s=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[-2]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[-1]))*180.0/np.pi),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
  # ax.set_title(key,fontsize=16)
  ax.legend(fontsize=10,ncol=3,loc="lower center")
  plt.savefig(plotFolder+"/energySpecTrigInclined"+str(suffix)+"scale"+str(energyScale)+"SelectCont"+str(containment)+".pdf",transparent=False,bbox_inches='tight')
  plt.close()

# plotInclinedEnergyFluxSelect(evtList,triggerListSelect,"log","fluxLog",0.0,containment=False)
plotInclinedEnergyFluxSelect(evtList,triggerListSelectDict,"log","fluxLog",0.0,containment=True)




def plotZenithFlux(eventList,triggerType,yscale,suffix,containment):
  """
  plots energy flux
  """
  if containment == True:
    # eventList = containedEvents(eventList,410)
    eventList = containedEvents(eventList,800)
  eventList = [ievt for ievt in eventList if abs(getattr(ievt,itrigger)-1)<0.01]
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(ncols=1,nrows=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  hitBins = np.linspace(14.0,17.0,31)
  totalRate = 0
  for nbin, binStart in enumerate(energyBinsShort[:-1]):
    lowEdge = energyBinsShort[nbin]
    highEdge = energyBinsShort[nbin+1]
    evtEnergyBin = [ievt for ievt in eventList if lowEdge <= ievt.energy < highEdge]
    zenithList = []
    neventList = []
    # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
    ncolor = next(colorIter)
    ax,histSum = plotZenithSteps(evtEnergyBin,ax,r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),ncolor=colorsCustom2[nbin])
    totalRate += histSum
  # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
  # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
  ax.set_xlabel(r"zenith [$^{\circ}$]", fontsize=22)
  ax.set_xticks(np.linspace(0,70,8))
  # ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20))
  ax.set_ylabel(r"log10(rate [Hz])", fontsize=20)
  # ax.set_yscale(yscale) 
  # ax.set_ylim(0,65)
  ax.xaxis.set_minor_locator(MultipleLocator(5))
  ax.set_xlim(0,65)
  ax.set_ylim(-4.5,0.5)
  # ax.set_xscale('log')
  ax.grid(True,alpha=0.7)
  ax.text(0.82,0.82,s=r"{0} rate:{1:.2f} Hz".format(triggerListSelectDict[triggerType],totalRate),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.text(0.9,0.9,s=r"{0:.1f} Hz".format(totalRate),horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, bbox=dict(boxstyle='round',facecolor='purple', alpha=0.1))
  # ax.set_title(key,fontsize=16)
  # ax.legend(fontsize=10,ncol=3,loc="lower center",columnspacing=0.4)
  ax.legend(fontsize=10,ncol=1,loc="upper left",columnspacing=0.4)
  plt.savefig(plotFolder+"/ZenithSpecTrig"+str(triggerListSelectDict[triggerType])+str(suffix)+"scaleCont"+str(containment)+".pdf",transparent=False,bbox_inches='tight')
  plt.close()

def plotZenithSteps(triggeredEvts,ax,legendLabel,ncolor):
  zenithBins = [np.arcsin(np.sqrt(izen)) for izen in sin2ZenBins]
  zenith = [ievt.zenith for ievt in triggeredEvts]
  weights = [ievt.H4aWeight for ievt in triggeredEvts]
  # weights_direct = [ievt.directWeight for ievt in triggeredEvts]
  hist,binEdge = np.histogram(zenith,zenithBins,weights=[w for w in weights])
  # binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
  print("sum of simulated rate",sum(hist))
  # if "total" in legendLabel:
    # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
  # ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
    # ax.set_ylim(10**-4,20)
  print("binedge",(binEdge[:-1]+binEdge[1:])/2.0)
  binCenter = (binEdge[:-1]+binEdge[1:])/2.0*180.0/np.pi
  H = [np.log10(h) for h in hist]   
  # ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
  ax.step(binCenter,H,"-",where="mid",lw=2.5,label=legendLabel+r", {:.4f} Hz".format(sum(hist)),color=ncolor,alpha=1)
  return ax,sum(hist)

# for itrigger in triggerListSelectDict.keys():
#   # plotZenithFlux(evtList,itrigger,"log","fluxLog",containment=False)
#   plotZenithFlux(evtList,itrigger,"log","fluxLog",containment=True)


# def plotEnergyFluxRatio(eventList1,eventList2,yscale,suffix,triggerType):
#   """
#   plots energy flux
#   """
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   hitBins = np.linspace(14.0,17.0,31)
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#     highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#     evtZenBin1 = [ievt for ievt in eventList1 if lowEdge <= ievt.zenith < highEdge]
#     evtZenBin2 = [ievt for ievt in eventList2 if lowEdge <= ievt.zenith < highEdge]
#     energyList = []
#     neventList = []
#     ncolor = next(colorIter)
#     # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#     ax = plotRatioSteps(evtZenBin1,evtZenBin2,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),ncolor=colorsCustom2[nbin])
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"$\rm \frac{{N_{{{0}}}}}{{N_{{{1}}}}}$".format(triggerType,"IceTopSMT"), fontsize=20)
#   ax.set_yscale(yscale) 
#   ax.set_ylim(10**(-0.5),10**1.8)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.6)
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3)
#   triggerType = re.sub(r"[^a-zA-Z0-9_ ]", "", triggerType)
#   plt.savefig(plotFolder+r"/energySpecRatioTrig"+str(triggerType)+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank6_5000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank6\_5000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank7_5000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank7\_5000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank8_5000-1)<0.01 ],
#   "log","fluxLog",triggerType=r"tank8\_5000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank9_5000-1)<0.01 ],
#   "log","fluxLog",triggerType=r"tank9\_5000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank10_5000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank10\_5000")

# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank6_4000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank6\_4000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank7_4000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank7\_4000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank8_4000-1)<0.01 ],
#   "log","fluxLog",triggerType=r"tank8\_4000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank9_4000-1)<0.01 ],
#   "log","fluxLog",triggerType=r"tank9\_4000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank10_4000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank10\_4000")

# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank6_3000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank6\_3000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank7_3000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank7\_3000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank8_3000-1)<0.01 ],
#   "log","fluxLog",triggerType=r"tank8\_3000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank9_3000-1)<0.01 ],
#   "log","fluxLog",triggerType=r"tank9\_3000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank10_3000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank10\_3000")

# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank6_2000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank6\_2000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank7_2000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank7\_2000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank8_2000-1)<0.01 ],
#   "log","fluxLog",triggerType=r"tank8\_2000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank9_2000-1)<0.01 ],
#   "log","fluxLog",triggerType=r"tank9\_2000")
# plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank10_2000-1)<0.01 ]
#   ,"log","fluxLog",triggerType=r"tank10\_2000")
# # plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank7Trig-1)<0.01 ],
# #   "log","fluxLog",triggerType="tank7")
# # plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.tank8Trig-1)<0.01 ],
# #   "log","fluxLog",triggerType="tank8")
# # plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.slc3Trig-1)<0.01 ]
# #   ,"log","fluxLog",triggerType="slc3")
# # plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.slc4Trig-1)<0.01 ],
# #   "log","fluxLog",triggerType="slc4")
# # plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.slc5Trig-1)<0.01 ],
# #   "log","fluxLog",triggerType="slc5")
# # plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.slc6Trig-1)<0.01 ]
# #   ,"log","fluxLog",triggerType="slc6")
# # plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.slc7Trig-1)<0.01 ],
# #   "log","fluxLog",triggerType="slc7")
# # plotEnergyFluxRatio([ievt for ievt in evtList if abs(ievt.HLC6_5000-1)<0.01],[ievt for ievt in evtList if abs(ievt.slc8Trig-1)<0.01 ],
# #   "log","fluxLog",triggerType="slc8")


# def plotInclinedEnergyFluxRatio(eventList,triggerList,yscale,suffix):
#   """
#   plots energy flux
#   """
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   hitBins = np.linspace(14.0,17.0,31)
#   lowEdge = np.arcsin(np.sqrt(sin2ZenBins[-2]))
#   highEdge = np.arcsin(np.sqrt(sin2ZenBins[-1]))
#   evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
#   evtZenBin1 = [ievt for ievt in evtZenBin if abs(getattr(ievt,"HLC6_5000")-1)<0.01]
#   for ntrig,itrigger in enumerate(triggerList[1:]):
#     evtZenBin2 = [ievt for ievt in evtZenBin if abs(getattr(ievt,itrigger)-1)<0.01]
#     energyList = []
#     neventList = []
#     ncolor = multiple_color[ntrig]
#     # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#     ax = plotRatioSteps(evtZenBin1,evtZenBin2,ax,str(itrigger),ncolor=ncolor)
#   # r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[-2]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[-1]))*180.0/np.pi),
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"$\rm \frac{{N_{{{0}}}}}{{N_{{{1}}}}}$".format(r"n\_tank","IceTopSMT"), fontsize=20)
#   ax.set_yscale(yscale) 
#   ax.set_ylim(10**(-0.5),10**1.8)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.6)
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,loc="upper left")
#   plt.savefig(plotFolder+r"/energySpecRatioTrigInclined"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# plotInclinedEnergyFluxRatio(evtList,triggerList,"log","fluxLog")

# def plotInclinedEnergyFluxRatio7(eventList,triggerList,yscale,suffix):
#   """
#   plots energy flux
#   """
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   hitBins = np.linspace(14.0,17.0,31)
#   lowEdge = np.arcsin(np.sqrt(sin2ZenBins[-2]))
#   highEdge = np.arcsin(np.sqrt(sin2ZenBins[-1]))
#   evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
#   evtZenBin1 = [ievt for ievt in evtZenBin if abs(getattr(ievt,"HLC6_5000")-1)<0.01]
#   for ntrig,itrigger in enumerate(triggerList[1:]):
#     evtZenBin2 = [ievt for ievt in evtZenBin if abs(getattr(ievt,itrigger)-1)<0.01]
#     energyList = []
#     neventList = []
#     ncolor = multiple_color[ntrig]
#     # ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#     ax = plotRatioSteps(evtZenBin1,evtZenBin2,ax,str(itrigger),ncolor=ncolor)
#   # r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[-2]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[-1]))*180.0/np.pi),
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"$\rm \frac{{N_{{{0}}}}}{{N_{{{1}}}}}$".format(r"n\_tank","IceTopSMT"), fontsize=20)
#   ax.set_yscale(yscale) 
#   ax.set_ylim(10**(-0.5),10**1.8)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.6)
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,loc="upper left")
#   plt.savefig(plotFolder+r"/energySpecRatioTrigInclined"+str(suffix)+"7.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# plotInclinedEnergyFluxRatio7(evtList,triggerList7,"log","fluxLog")


def plotTrigEfficiency(evtList,energyBins,triggerType,containment):
  '''
  plots trigger efficiency in different zenith bins
  '''
  print("plotting trigger efficiency for ",triggerType)
  if containment == True:
    # evtList = containedEvents(evtList,640)
    evtList = containedEvents(evtList,410)
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # colorIter = iter(colorsList)
  # energyBins = 10**np.linspace(5.0, 8.0, 31)
  energyBins = 10**np.linspace(5.0, 7.8, 29)
  for nbin, binStart in enumerate(sin2ZenBins[:-1]):
    lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
    highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
    evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
    energyList = []
    efficiencyList = []
    ncolor = colorsCustom2[nbin]
    for ebin, ebinStart in enumerate(energyBins[:-1]):
      lowEdge_E = energyBins[ebin]
      highEdge_E = energyBins[ebin+1]
      evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
      # totalEvts = len(evtEBin)
      weights = [ievt.H4aWeight for ievt in evtEBin]
      totalEvts = len(evtEBin)
      # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
      # triggerList = [ievt for ievt in evtEBin if abs(getattr(ievt,triggerType)-1)<0.01]
      triggerList = [ievt for ievt in evtEBin if abs(getattr(ievt,triggerType)-0)>0.01]
      trigEff = triggerEfficiency(len(triggerList),totalEvts)
      efficiencyList.append(trigEff)
      energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))   
    ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  ax.set_ylabel(r"trigger efficiency", fontsize=22)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  ax.text(0.78,0.1,s=r"trig:{0}".format(triggerListSelectDict[triggerType]),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  ax.text(0.78,0.2,s=r"snow:{0}".format("2021/03"),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_xscale('log')
  ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  ax.set_ylim(0,1.05)
  # ax.set_xlim(14,17)
  ax.set_xlim(14,16.8)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  ax.grid(True,alpha=0.5)
  # l1=ax.legend(loc="upper left",fontsize=12)
  l1=ax.legend(loc="upper left",fontsize=12)
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"0.98")
  # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  # ax.add_artist(l1)
  # ax.add_artist(l2)
  plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"Efficiency.pdf",transparent=False,bbox_inches='tight')
  plt.close()


# plotTrigEfficiency(evtList,energyBins,triggerType="tank1",containment=True)
# plotTrigEfficiency(evtList,energyBins,triggerType="sta1",containment=True)
###########################################################
# for itrigger in triggerList:
#   plotTrigEfficiency(evtList,energyBins,triggerType=itrigger,containment=True)
#########################################################
for itrigger in triggerListSelectDict.keys():
  plotTrigEfficiency(evtList,energyBins,triggerType=itrigger,containment=True)
# ########################################################
# #########################################################


# def plotInclinedTrigEfficiency(evtList,energyBins,triggerTypes,containment):
#   '''
#   plots trigger efficiency in different zenith bins
#   '''
#   print("plotting trigger efficiency for ",triggerTypes)
#   if containment == True:
#     # evtList = containedEvents(evtList,640)
#     evtList = containedEvents(evtList,410)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   # colorIter = iter(colorsList)
#   for ntrig,itrigger in enumerate(triggerTypes):
#       lowEdge = np.arcsin(np.sqrt(sin2ZenBins[-2]))
#       highEdge = np.arcsin(np.sqrt(sin2ZenBins[-1]))
#       evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
#       energyList = []
#       efficiencyList = []
#       # ncolor = colorsCustom2[ntrig]
#       ncolor = multiple_color[ntrig]
#       for ebin, ebinStart in enumerate(energyBins[:-1]):
#         lowEdge_E = energyBins[ebin]
#         highEdge_E = energyBins[ebin+1]
#         evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
#         # totalEvts = len(evtEBin)
#         weights = [ievt.H4aWeight for ievt in evtEBin]
#         totalEvts = len(evtEBin)
#         # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
#         triggerList = [ievt for ievt in evtEBin if abs(int(getattr(ievt,str(itrigger)))-1)<0.01]
#         trigEff = triggerEfficiency(len(triggerList),totalEvts)
#         efficiencyList.append(trigEff)
#         energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))   
#       ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=str(itrigger),alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"trigger efficiency", fontsize=22)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   ax.text(0.78,0.1,s=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[-2]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[-1]))*180.0/np.pi),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#   # ax.set_xscale('log')
#   ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
#   ax.set_ylim(0,1.05)
#   ax.set_xlim(14,17)
#   # ax.yaxis.set_minor_locator(MultipleLocator(100))
#   ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#   ax.grid(True,alpha=0.5)
#   l1=ax.legend(loc="upper left",fontsize=12)
#   point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='gray', marker='',markersize=5, label=r"0.98")
#   # l2 = ax.legend(handles=[point_dash],loc="center left",fontsize=13,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
#   l2 = ax.legend(handles=[point_dash],fontsize=13,framealpha=0.1,handlelength=1.4,handletextpad=0.5,bbox_to_anchor=(0.85,0.25),bbox_transform=ax.transAxes,prop={"family":"serif","size":13})
#   ax.add_artist(l1)
#   ax.add_artist(l2)
#   plt.savefig(plotFolder+"/trigInclinedcont"+str(containment)+"Efficiency.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# plotInclinedTrigEfficiency(evtList,energyBins,triggerTypes=triggerList,containment=True)

# def plotInclinedTrigEfficiency7(evtList,energyBins,triggerTypes,containment):
#   '''
#   plots trigger efficiency in different zenith bins
#   '''
#   print("plotting trigger efficiency for ",triggerTypes)
#   if containment == True:
#     # evtList = containedEvents(evtList,640)
#     evtList = containedEvents(evtList,410)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsCustom+colorsCustom)
#   # colorIter = iter(colorsList)
#   for ntrig,itrigger in enumerate(triggerTypes):
#       lowEdge = np.arcsin(np.sqrt(sin2ZenBins[-2]))
#       highEdge = np.arcsin(np.sqrt(sin2ZenBins[-1]))
#       evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
#       energyList = []
#       efficiencyList = []
#       # ncolor = colorsCustom2[ntrig]
#       ncolor = multiple_color[ntrig]
#       for ebin, ebinStart in enumerate(energyBins[:-1]):
#         lowEdge_E = energyBins[ebin]
#         highEdge_E = energyBins[ebin+1]
#         evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
#         # totalEvts = len(evtEBin)
#         weights = [ievt.H4aWeight for ievt in evtEBin]
#         totalEvts = len(evtEBin)
#         # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
#         triggerList = [ievt for ievt in evtEBin if abs(int(getattr(ievt,str(itrigger)))-1)<0.01]
#         trigEff = triggerEfficiency(len(triggerList),totalEvts)
#         efficiencyList.append(trigEff)
#         energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))   
#       ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=str(itrigger),alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"trigger efficiency", fontsize=22)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   ax.text(0.78,0.1,s=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[-2]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[-1]))*180.0/np.pi),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#   # ax.set_xscale('log')
#   ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
#   ax.set_ylim(0,1.05)
#   ax.set_xlim(14,17)
#   # ax.yaxis.set_minor_locator(MultipleLocator(100))
#   ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#   ax.grid(True,alpha=0.5)
#   l1=ax.legend(loc="upper left",fontsize=12)
#   point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='gray', marker='',markersize=5, label=r"0.98")
#   # l2 = ax.legend(handles=[point_dash],loc="center left",fontsize=13,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
#   l2 = ax.legend(handles=[point_dash],fontsize=13,framealpha=0.1,handlelength=1.4,handletextpad=0.5,bbox_to_anchor=(0.85,0.25),bbox_transform=ax.transAxes,prop={"family":"serif","size":13})
#   ax.add_artist(l1)
#   ax.add_artist(l2)
#   plt.savefig(plotFolder+"/trigInclinedcont"+str(containment)+"Efficiency7.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# plotInclinedTrigEfficiency7(evtList,energyBins,triggerTypes=triggerList7,containment=True)

def plotVerticalTrigEfficiencySelect(evtList,energyBins,triggerTypes,containment,wilson):
  '''
  plots trigger efficiency in different zenith bins
  '''
  print("plotting trigger efficiency for ",triggerTypes)
  if containment == True:
    # evtList = containedEvents(evtList,640)
    evtList = containedEvents(evtList,410)
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # colorIter = iter(colorsList)
  for ntrig,itrigger in enumerate(triggerTypes):
      lowEdge = np.arcsin(np.sqrt(sin2ZenBins[0]))
      highEdge = np.arcsin(np.sqrt(sin2ZenBins[1]))
      evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
      energyList = []
      efficiencyList = []
      wilsonMeanList = []
      errorListLow = []
      errorListHigh = []
      # ncolor = colorsCustom2[ntrig]
      ncolor = multiple_color[ntrig]
      for ebin, ebinStart in enumerate(energyBins[:-1]):
        lowEdge_E = energyBins[ebin]
        highEdge_E = energyBins[ebin+1]
        evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
        # totalEvts = len(evtEBin)
        weights = [ievt.H4aWeight for ievt in evtEBin]
        totalEvts = len(evtEBin)
        # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
        triggerList = [ievt for ievt in evtEBin if abs(int(getattr(ievt,str(itrigger)))-1)<0.01]
        trigEff = triggerEfficiency(len(triggerList),totalEvts)
        efficiencyList.append(trigEff)
        #####################################################################################
        #binomial interval
        # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68)
        #####################################################################################
        wilsonM = WilsonMean(nPass=len(triggerList), nFail=totalEvts-len(triggerList))
        wilsonErr = WilsonError(nPass=len(triggerList), nFail=totalEvts-len(triggerList))
        ErrorH = wilsonM + wilsonErr - trigEff
        ErrorL = trigEff - (wilsonM - wilsonErr)
        # print("wilson calculation",totalEvts,len(triggerList),trigEff,wilsonErr)
        wilsonMeanList.append(wilsonM)
        errorListLow.append(ErrorL)
        errorListHigh.append(ErrorH)
        energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))
      # ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",ls="-",lw = 2.5,c=ncolor,label=str(triggerListSelectDict[itrigger]),alpha=1)
      # ax.errorbar(energyList,wilsonMeanList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="none",lw=2.5,c=ncolor,label=str(triggerListSelectDict[itrigger]),alpha=1)
      ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",ls="-",lw=2.5,c=ncolor,label=str(triggerListSelectDict[itrigger]),alpha=1)
      ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  ax.set_ylabel(r"trigger efficiency", fontsize=22)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  ax.text(0.78,0.1,s=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[0]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[1]))*180.0/np.pi),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  ax.text(0.12,0.7,s=r"snow:{0}".format("2021/03"),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_xscale('log')
  ax.text(0.75,0.5,s="IceCube Preliminary",color="red",size=14,fontWeight='bold',horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  ax.set_ylim(0.0,1.05)
  # ax.set_ylim(0.9,1.01)
  # ax.set_xlim(14,17)
  ax.set_xlim(14,16.8)
  # ax.set_xlim(15.8,17)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  ax.grid(True,alpha=0.5)
  l1=ax.legend(loc="upper left",fontsize=12)
  # l1=ax.legend(loc="upper left",fontsize=12)
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='gray', marker='',markersize=5, label=r"0.98")
  # l2 = ax.legend(handles=[point_dash],loc="center left",fontsize=13,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  l2 = ax.legend(handles=[point_dash],fontsize=13,framealpha=0.1,handlelength=1.4,handletextpad=0.5,bbox_to_anchor=(0.85,0.25),bbox_transform=ax.transAxes,prop={"family":"serif","size":13})
  ax.add_artist(l1)
  ax.add_artist(l2)
  plt.savefig(plotFolder+"/trigVertcont"+str(containment)+"EfficiencySelect.pdf",transparent=False,bbox_inches='tight')
  plt.close()

# plotVerticalTrigEfficiencySelect(evtList,energyBins,triggerTypes=triggerListSelect,containment=True)
# plotVerticalTrigEfficiencySelect(evtList,10**np.linspace(5.0, 8.0, 30),triggerTypes=triggerListSelectDict.keys(),containment=True)
plotVerticalTrigEfficiencySelect(evtList,10**np.linspace(5.0, 7.8, 29),triggerTypes=triggerListSelectDict.keys(),containment=True,wilson=True)


def plotInclinedTrigEfficiencySelect(evtList,energyBins,triggerTypes,containment,wilson):
  '''
  plots trigger efficiency in different zenith bins
  '''
  print("plotting trigger efficiency for ",triggerTypes)
  if containment == True:
    # evtList = containedEvents(evtList,640)
    evtList = containedEvents(evtList,410)
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # colorIter = iter(colorsList)
  for ntrig,itrigger in enumerate(triggerTypes):
      lowEdge = np.arcsin(np.sqrt(sin2ZenBins[-2]))
      highEdge = np.arcsin(np.sqrt(sin2ZenBins[-1]))
      evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
      energyList = []
      efficiencyList = []
      wilsonMeanList = []
      errorListLow = []
      errorListHigh = []
      # ncolor = colorsCustom2[ntrig]
      ncolor = multiple_color[ntrig]
      for ebin, ebinStart in enumerate(energyBins[:-1]):
        lowEdge_E = energyBins[ebin]
        highEdge_E = energyBins[ebin+1]
        evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
        # totalEvts = len(evtEBin)
        weights = [ievt.H4aWeight for ievt in evtEBin]
        totalEvts = len(evtEBin)
        # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
        triggerList = [ievt for ievt in evtEBin if abs(int(getattr(ievt,str(itrigger)))-1)<0.01]
        trigEff = triggerEfficiency(len(triggerList),totalEvts)
        efficiencyList.append(trigEff)
        #binomial interval
        # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68)
        #####################################################################################
        wilsonM = WilsonMean(nPass=len(triggerList), nFail=totalEvts-len(triggerList))
        wilsonErr = WilsonError(nPass=len(triggerList), nFail=totalEvts-len(triggerList))
        ErrorH = wilsonM + wilsonErr - trigEff
        ErrorL = trigEff - (wilsonM - wilsonErr)
        # print("wilson calculation",totalEvts,len(triggerList),trigEff,wilsonErr)
        wilsonMeanList.append(wilsonM)
        errorListLow.append(ErrorL)
        errorListHigh.append(ErrorH)
        # print("efficiency bug",wilsonM, wilsonErr)
        energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))
      # ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",ls="-",lw = 2.5,c=ncolor,label=str(triggerListSelectDict[itrigger]),alpha=1)
      ax.errorbar(energyList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",ls="-",lw=2.5,label=str(triggerListSelectDict[itrigger]),c=ncolor,alpha=1)
      ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
  ax.set_ylabel(r"trigger efficiency", fontsize=22)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  ax.text(0.78,0.1,s=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[-2]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[-1]))*180.0/np.pi),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  ax.text(0.15,0.3,s=r"snow:{0}".format("2021/03"),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  ax.text(0.20,0.5,s="IceCube Preliminary",color="red",size=14,fontWeight='bold',horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_xscale('log')
  ax.axhline(y=0.98,xmin=0,xmax=1,color="gray",linestyle="--",lw=2.0)
  ax.set_ylim(0.0,1.05)
  # ax.set_ylim(0.9,1.01)
  # ax.set_xlim(14,17)
  ax.set_xlim(14,16.8)
  # ax.set_xlim(15.8,17)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  ax.xaxis.set_minor_locator(MultipleLocator(0.1))
  ax.grid(True,alpha=0.5)
  l1=ax.legend(loc="lower left",fontsize=12)
  # l1=ax.legend(loc="upper left",fontsize=12)
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='gray', marker='',markersize=5, label=r"0.98")
  # l2 = ax.legend(handles=[point_dash],loc="center left",fontsize=13,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  l2 = ax.legend(handles=[point_dash],fontsize=13,framealpha=0.1,handlelength=1.4,handletextpad=0.5,bbox_to_anchor=(0.85,0.25),bbox_transform=ax.transAxes,prop={"family":"serif","size":13})
  ax.add_artist(l1)
  ax.add_artist(l2)
  plt.savefig(plotFolder+"/trigInclinedcont"+str(containment)+"EfficiencySelect.pdf",transparent=False,bbox_inches='tight')
  plt.close()

# plotInclinedTrigEfficiencySelect(evtList,energyBins,triggerTypes=triggerListSelect,containment=True)
# plotInclinedTrigEfficiencySelect(evtList,10**np.linspace(5.0, 8.0, 30),triggerTypes=triggerListSelectDict.keys(),containment=True)
plotInclinedTrigEfficiencySelect(evtList,10**np.linspace(5.0, 7.8, 29),triggerTypes=triggerListSelectDict.keys(),containment=True,wilson=True)


def plotTrigEfficiencyZenith(evtList,energyBins,triggerType,containment,wilson):
  '''
  plots trigger efficiency in different zenith bins
  '''
  print("plotting trigger efficiency for ",triggerType)
  if containment == True:
    # evtList = containedEvents(evtList,640)
    evtList = containedEvents(evtList,410)
  # sin2ZenBins = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.822]
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  colorIter = iter(colorsCustom+colorsCustom)
  # colorIter = iter(colorsList)
  # zorders = [0,1,2,3,4,5,6]
  zorders = [60,50,40,30,20,10,9]
  for nbin, binStart in enumerate(energyBinsShort[:-1]):
    lowEdge = energyBinsShort[nbin]
    highEdge = energyBinsShort[nbin+1]
    evtEnergyBin = [ievt for ievt in evtList if lowEdge <= ievt.energy < highEdge]
    zenithList = []
    efficiencyList = []
    wilsonMeanList = []
    errorListLow = []
    errorListHigh = []
    ncolor = colorsCustom2[nbin]
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in evtEnergyBin if lowEdge_Z <= np.sin(ievt.zenith)**2 < highEdge_Z]
      weights = [ievt.H4aWeight for ievt in evtZBin]
      totalEvts = len(evtZBin)
      # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
      triggerList = [ievt for ievt in evtZBin if abs(getattr(ievt,triggerType)-1)<0.01]
      trigEff = triggerEfficiency(len(triggerList),totalEvts)
      #########################################################
      # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68) #binomial interval
      ##########################################################
      #binomial interval
      # wilsonErr = binomial_proportion(nsel=len(triggerList), ntot=totalEvts,coverage=0.68)
      #####################################################################################
      wilsonM = WilsonMean(nPass=len(triggerList), nFail=totalEvts-len(triggerList))
      wilsonErr = WilsonError(nPass=len(triggerList), nFail=totalEvts-len(triggerList))
      ErrorH = wilsonM + wilsonErr - trigEff
      ErrorL = trigEff - (wilsonM - wilsonErr)
      # print("wilson calculation",totalEvts,len(triggerList),trigEff,wilsonErr)
      wilsonMeanList.append(wilsonM)
      errorListLow.append(ErrorL)
      errorListHigh.append(ErrorH)
      efficiencyList.append(trigEff)
      zenithList.append((np.arcsin(np.sqrt(lowEdge_Z))+np.arcsin(np.sqrt(highEdge_Z)))/2.0*180/np.pi)
    # ax.plot(zenithList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
    # ax.errorbar(zenithList,wilsonMeanList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
    # ax.errorbar(zenithList,wilsonMeanList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="none",ls = 2.5,c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
    ax.errorbar(zenithList,efficiencyList,yerr=np.asarray([errorListLow,errorListHigh]),fmt="o",ls="-",lw = 2.5,c=ncolor,label=r"{0:.1f}-{1:.1f}".format(np.log10(energyBinsShort[nbin])+9,np.log10(energyBinsShort[nbin+1])+9),zorder=zorders[nbin],alpha=1)
    # ax.plot(zenithList,efficiencyList,".",ls='-',lw = 2.5,c=ncolor,zorder=zorders[nbin],alpha=1)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
  ax.set_xlabel(r"zenith [$^{\circ}$]", fontsize=22)
  ax.set_ylabel(r"trigger efficiency", fontsize=22)
  ax.text(0.20,0.25,s="IceCube Preliminary",color="red",size=14,fontWeight='bold',horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
  # ax.text(0.78,0.1,s=r"trig:{0}".format(triggerListSelectDict[triggerType]),size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.set_xscale('log')
  ax.set_xticks(np.linspace(0,70,8))
  yline=0.98
  ax.axhline(y=yline,xmin=0,xmax=1,color="gray",linestyle="--",zorder=1,lw=2.0)
  ax.set_ylim(0,1.05)
  ax.set_xlim(0,65)
  # ax.yaxis.set_minor_locator(MultipleLocator(100))
  ax.xaxis.set_minor_locator(MultipleLocator(5))
  ax.grid(True,alpha=0.5)
  l1=ax.legend(loc="upper left",title=r"log$_{10}$ [E/eV]",fontsize=16).set_zorder(10000)
  # l1=ax.legend(loc="upper left",title=r"log$_{10}$ [E/eV]",title_fontsize=14,fontsize=14).set_zorder(10000)
  # l1=ax.legend(loc="upper left",fontsize=14).set_zorder(10000)
  # plt.setp(l1.get_title(),fontsize=14)
  # plt.legend.set_title(r"log$_{10}$ (E/eV)",prop={'size':14})
  point_dash = mlines.Line2D([], [], linestyle='--',lw=2.0,color='black', marker='',markersize=5, label=r"{yline:.2f}")
  # l2 = ax.legend(handles=[point_dash],loc="upper left",fontsize=20,framealpha=0.1,handlelength=1.4,handletextpad=0.5)
  # ax.add_artist(l1)
  # ax.add_artist(l2)
  plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"EfficiencyZenith.pdf",transparent=False,bbox_inches='tight')
  plt.close()


for itrigger in triggerListSelectDict.keys():
  # plotTrigEfficiencyZenith(evtList,energyBins,triggerType=itrigger,containment=True,wilson=False)
  plotTrigEfficiencyZenith(evtList,energyBins,triggerType=itrigger,containment=True,wilson=True)



# def delta_t_hist_zen_bins(evtList,deltaType):
#   '''
#   plots histogram of delataT
#   '''
#   print("plotting deltaT histogram",deltaType)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   # colorIter = iter(colorsCustom)
#   colorIter = iter(colorsList)
#   histBins = np.linspace(0,6000,6000)
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#     highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#     evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
#     if deltaType == "slc3":
#       deltaTList = [ievt.deltaT3SLC for ievt in evtZenBin]
#     elif deltaType == "slc4":
#       deltaTList = [ievt.deltaT4SLC for ievt in evtZenBin]
#     elif deltaType == "slc5":
#       deltaTList = [ievt.deltaT5SLC for ievt in evtZenBin]
#     elif deltaType == "slc6":
#       deltaTList = [ievt.deltaT6SLC for ievt in evtZenBin]
#     elif deltaType == "slc7":
#       deltaTList = [ievt.deltaT7SLC for ievt in evtZenBin]
#     elif deltaType == "slc8":
#       deltaTList = [ievt.deltaT8SLC for ievt in evtZenBin]
#     elif deltaType == "tank3":
#       deltaTList = [ievt.deltaT3Tank for ievt in evtZenBin]
#     elif deltaType == "tank4":
#       deltaTList = [ievt.deltaT4Tank for ievt in evtZenBin]
#     elif deltaType == "tank5":
#       deltaTList = [ievt.deltaT5Tank for ievt in evtZenBin]
#     elif deltaType == "tank6":
#       deltaTList = [ievt.deltaT6Tank for ievt in evtZenBin]
#     elif deltaType == "tank7":
#       deltaTList = [ievt.deltaT7Tank for ievt in evtZenBin]
#     elif deltaType == "tank8":
#       deltaTList = [ievt.deltaT8Tank for ievt in evtZenBin]
#     print("deltaTList",len(deltaTList))
#     ax.hist(deltaTList,bins=histBins,histtype="step",label=r"$\theta$ = {0:.0f}$^{{\circ}}$-{1:.0f}$^{{\circ}}$, {2:d} evts".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi,len(deltaTList)),lw=2.5,alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"$\Delta_t$ [ns]", fontsize=20)
#   ax.set_ylabel(r"count", fontsize=20)
#   ax.set_yscale('log')
#   ax.set_xlim(0,6000)
#   ax.set_ylim(None,10**3.1)
#   ax.grid(True,alpha=0.7)
#   # ax.yaxis.set_minor_locator(MultipleLocator(100))
#   # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#   ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/deltaT"+str(deltaType)+"ZenBinsdelta_t.pdf",transparent=False,bbox_inches='tight')
#   plt.close()


# delta_t_hist_zen_bins(evtList,"tank3")
# delta_t_hist_zen_bins(evtList,"tank4")
# delta_t_hist_zen_bins(evtList,"tank5")
# delta_t_hist_zen_bins(evtList,"tank6")
# delta_t_hist_zen_bins(evtList,"tank7")
# delta_t_hist_zen_bins(evtList,"tank8")
# delta_t_hist_zen_bins(evtList,"slc3")
# delta_t_hist_zen_bins(evtList,"slc4")
# delta_t_hist_zen_bins(evtList,"slc5")
# delta_t_hist_zen_bins(evtList,"slc6")
# delta_t_hist_zen_bins(evtList,"slc7")
# delta_t_hist_zen_bins(evtList,"slc8")


##################################################################################################################################
##################################################################################################################################


# def addHistogram(x,ax,suffix):
#     xbins = np.linspace(min(x),max(x),200)
#     ax.hist(x,bins=xbins,histtype="step",lw=2.5,alpha=1)
#     return ax

# def plotDistanceHist(evtList,triggerType):
#   '''
#   plots core distance of triggered shower
#   '''
#   energyBins = 10**(np.linspace(5.0,8.0,7))
#   evtList = selectTriggered(evtList,triggerType)
#   for ebin, ebinStart in enumerate(energyBins[:-1]):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     colorIter = iter(colorsList)
#     lowEdge_E = energyBins[ebin]
#     highEdge_E = energyBins[ebin+1]
#     evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
#     ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((energyBins[ebin])*10**9),np.log10((energyBins[ebin+1])*10**9)))
#     for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#       # colorIter = iter(colorsCustom)
#       lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#       highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#       evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
#       distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenBin]
#       if len(distanceList)>2:
#         xbins = np.linspace(min(distanceList),max(distanceList),200)
#       else:
#         xbins = np.linspace(0,1000,200)
#       ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
#         r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),lw=2.5,alpha=1)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#     ax.set_xlabel(r"core distance [m]", fontsize=22)
#     ax.set_ylabel(r"count", fontsize=22)
#     # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#     # ax.set_xscale('log')
#     # ax.set_ylim(0,1.05)
#     # ax.set_xlim(14,17)
#     # ax.yaxis.set_minor_locator(MultipleLocator(100))
#     # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#     ax.grid(True,alpha=0.6)
#     ax.legend(fontsize=12)
#     plt.savefig(plotFolder+"/distanceHistTrig"+str(triggerType)+"energy{0:.1f}.pdf".format(np.log10(energyBins[ebin])+9),transparent=False,bbox_inches='tight')
#     plt.close()
# plotDistanceHist(evtList,triggerType="sta1")
# plotDistanceHist(evtList,triggerType="slc3")
# plotDistanceHist(evtList,triggerType="sta3")

# def plotDistanceHistGiven(evtList,triggerType,energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[400,500]):
#   '''
#   plots core distance of triggered shower
#   '''
#   energyBins = 10**(np.linspace(5.0,8.0,7))
#   # evtList = selectTriggered(evtList,triggerType)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   colorIter = iter(colorsList)
#   lowEdge_E = energyLim[0]
#   highEdge_E = energyLim[1]
#   evtEBin = [ievt for ievt in evtList if lowEdge_E <= ievt.energy < highEdge_E]
#   ax.set_title(r"log10(E[eV]):{0:.1f}-{1:.1f}".format(np.log10((lowEdge_E)*10**9),np.log10((highEdge_E)*10**9)))
#   lowEdge = zenLim[0]
#   highEdge = zenLim[1]
#   evtZenBin = [ievt for ievt in evtEBin if lowEdge <= ievt.zenith < highEdge]
#   evtZenDistBin = [ievt for ievt in evtZenBin if distanceLim[0] <= np.sqrt(ievt.coreX**2+ievt.coreY**2) < distanceLim[1]]
#   distanceList = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenDistBin]
#   xbins = np.linspace(min(distanceList),max(distanceList),30)
#   evtZenDistBinTriggerd = selectTriggered(evtZenDistBin,triggerType)
#   distanceListTrig = [np.sqrt(ievt.coreX**2+ievt.coreY**2) for ievt in evtZenDistBinTriggerd]
#   ax.hist(distanceList,bins=xbins,histtype="step",color=next(colorIter),label=
#     r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(zenLim[0]*180.0/np.pi,zenLim[1]*180.0/np.pi),lw=2.5,alpha=0.4)
#   ax.hist(distanceListTrig,bins=xbins,histtype="step",color=next(colorIter),label=
#     r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(zenLim[0]*180.0/np.pi,zenLim[1]*180.0/np.pi)+" "+str(triggerType),lw=2.5,alpha=0.4)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"core distance [m]", fontsize=22)
#   ax.set_ylabel(r"count", fontsize=22)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   # ax.set_xscale('log')
#   # ax.set_ylim(0,1.05)
#   # ax.set_xlim(14,17)
#   # ax.yaxis.set_minor_locator(MultipleLocator(100))
#   # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#   ax.grid(True,alpha=0.6)
#   ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/distanceHistTrig"+str(triggerType)+"GivenEnergy{0:.1f}.pdf".format(lowEdge_E+9),transparent=False,bbox_inches='tight')
#   plt.close()

# plotDistanceHistGiven(evtList,triggerType="sta3",energyLim=[10**6.0,10**6.5],zenLim=[np.arcsin(np.sqrt(0.0)),np.arcsin(np.sqrt(0.1))],distanceLim=[442,446])


# # plotCoreScatter(hdf5NullList)

# def plotZenithHist(hdfFileList):
#   zenList,ptypeList,energyList = getZenithTypeEnergy(hdfFileList)
#   # zenList = [np.rad2deg(i) for i in zenList]
#   zenList = [np.sin(i)**2 for i in zenList]
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   ax.hist(zenList,label="radius",alpha=1)
#   # ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
#   ax.set_xlabel(r"log(Energy[GeV])", fontsize=20)
#   ax.set_ylabel(r"disc radius [m]", fontsize=20)
#   ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
#   # ax.set_yscale('log')
#   ax.grid(True,alpha=0.2)
#   ax.set_ylim(0,None)
#   # ax.legend(fontsize=14)
#   # ax.legend(fontsize=14,ncol=2)
#   plt.savefig(plotFolder+"zenHist.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # plotZenithHist(hdf5NullList)


# def plotEnergyFlux2(eventList,yscale,suffix,energyScale,triggerType):
#   """
#   plots energy flux
#   """
#   # triggeredEvts = [ievt.ITSTA3_filter for ievt in eventList]
#   # print("triggers",triggeredEvts)
#   if triggerType == "sta3":
#     triggeredEvts =[ievt for ievt in eventList if abs(ievt.ITSMTTriggered-1)<0.01]
#   elif triggerType == "sta1":
#     triggeredEvts =[ievt for ievt in eventList if abs(ievt.STA1Trigger-1)<0.01]
#   elif triggerType == "slc3":
#     triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc3Trig-1)<0.01]
#   elif triggerType == "slc4":
#     triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc4Trig-1)<0.01]
#   elif triggerType == "slc5":
#     triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc5Trig-1)<0.01]
#   energy = [ievt.energy for ievt in triggeredEvts]
#   weights = [ievt.H4aWeight for ievt in triggeredEvts]
#   weights_direct = [ievt.directWeight for ievt in triggeredEvts]
#   # print("trigger",len(triggeredEvts))
#   print("energy",energy[:5])
#   print("weights",weights[:5])
#   print("direct",weights_direct[:5])
#   energy = np.log10(energy)+9.0
#   # print("lengths",len(weights),weights,len(energy))
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   hitBins = np.linspace(14.0,17.0,31)
#   hist,binEdge = np.histogram(energy,hitBins,weights=[w for w in weights])
#   # binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
#   print("sum of simulated rate",sum(hist))
#   ax.text(0.05, 0.95, r"sim. rate:{:.1f} Hz".format(sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
#   binCenter = (binEdge[:-1]+binEdge[1:])/2.0
#   if str(energyScale) == "0.0":
#     H = [h for h in hist]
#     ax.set_ylim(10**-4,20)
#   else:
#     H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
#   # ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
#   ax.step(binCenter,H,"-",where="mid",lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20)
#   ax.set_yscale(yscale) 
#   # ax.set_ylim(None,10**5)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.2)
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()
# # plotEnergyFlux(energy,weights,adjustedWeights,"log","flux")
# # plotEnergyFlux(evtList,"linear","fluxlinear",1.8,triggerType="sta3")
# # plotEnergyFlux2(evtList,"linear","fluxlinear",0.0,triggerType="sta3")
# # plotEnergyFlux2(evtList,"linear","fluxlinear",0.0,triggerType="sta1")
# # plotEnergyFlux2(evtList,"log","fluxLog",0.0,triggerType="sta1")
# # plotEnergyFlux2(evtList,"log","fluxLog",0.0,triggerType="sta3")
# # plotEnergyFlux2(evtList,"log","fluxLog",1.8,triggerType="sta1")
# # plotEnergyFlux2(evtList,"log","fluxLog",1.8,triggerType="sta3")




# def plotEnergyFlux(eventList,yscale,suffix,energyScale,triggerType,plotTotal=True,sepZenBins=True):
#   """
#   plots energy flux
#   """
#   # triggeredEvts = [ievt.ITSTA3_filter for ievt in eventList]
#   # print("triggers",triggeredEvts)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   hitBins = np.linspace(14.0,17.0,31)
#   if sepZenBins == False:
#     if triggerType == "sta3":
#       triggeredEvts =[ievt for ievt in eventList if abs(ievt.ITSMTTriggered-1)<0.01]
#     elif triggerType == "sta1":
#       triggeredEvts =[ievt for ievt in eventList if abs(ievt.STA1Trigger-1)<0.01]
#     elif triggerType == "slc3":
#       triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc3Trig-1)<0.01]
#     elif triggerType == "slc4":
#       triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc4Trig-1)<0.01]
#     elif triggerType == "slc5":
#       triggeredEvts =[ievt for ievt in eventList if abs(ievt.slc5Trig-1)<0.01]
#     print("debug no of triggered events",len(triggeredEvts),triggerType)
#     if plotTotal==True:
#       ax = plotSteps(eventList,ax,"{} events".format("total"),energyScale)
#     ax = plotSteps(triggeredEvts,ax,"{} trigger".format(triggerType),energyScale)
#   elif sepZenBins == True:
#     for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#       lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#       highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#       evtZenBin = [ievt for ievt in eventList if lowEdge <= ievt.zenith < highEdge]
#       energyList = []
#       neventList = []
#       if triggerType == "sta3":
#         triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.ITSMTTriggered-1)<0.01]
#       elif triggerType == "sta1":
#         triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.STA1Trigger-1)<0.01]
#       elif triggerType == "slc3":
#         triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.slc3Trig-1)<0.01]
#       elif triggerType == "slc4":
#         triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.slc4Trig-1)<0.01]
#       elif triggerType == "slc5":
#         triggeredEvts =[ievt for ievt in evtZenBin if abs(ievt.slc5Trig-1)<0.01]
#       else:
#         print("unkbown trigger type",triggerType)
#       if plotTotal==True:
#         ax = plotSteps(evtZenBin,ax,"{}".format("total"),energyScale)
#       ax = plotSteps(triggeredEvts,ax,r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),energyScale)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"$E^{{{0:.1f}}}$ rate [Hz]".format(energyScale), fontsize=20)
#   ax.set_yscale(yscale) 
#   # ax.set_ylim(None,10**5)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.2)
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/energySpecTrig"+str(triggerType)+str(suffix)+"scale"+str(energyScale)+"Total"+str(plotTotal)+"ZenBin"+str(sepZenBins)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()
# def plotEffectiveArea(evtList,energyBins,weighting):
#   '''
#   plots trigger efficiency in different zenith bins
#   '''
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   # colorIter = iter(colorsCustom)
#   colorIter = iter(colorsList)
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#     highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#     evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
#     energyList = []
#     sta3_efficiencyList = []
#     sta1_efficiencyList = []
#     for ebin, ebinStart in enumerate(energyBins[:-1]):
#       lowEdge_E = energyBins[ebin]
#       highEdge_E = energyBins[ebin+1]
#       evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
#       # totalEvts = len(evtEBin)
#       weights = [ievt.H4aWeight for ievt in evtEBin]
#       # totalEvts = np.sum(weights)
#       totalEvts = len(evtEBin)
#       # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
#       sta3 = [ievt.ITSMTTriggered for ievt in evtEBin]
#       # sta3 = [ievt.ITSTA5_filter for ievt in evtEBin]
#       sta1 = [ievt.STA1Trigger for ievt in evtEBin]
#       sta3_trigEff = effectiveArea(sum(sta3),totalEvts,showerArea(np.log10(lowEdge_E)))
#       sta1_trigEff = triggerEfficiency(sum(sta1),totalEvts)
#       sta3_efficiencyList.append(sta3_trigEff)
#       sta1_efficiencyList.append(sta1_trigEff)
#       energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9)) 
#     ax.plot(energyList,sta3_efficiencyList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"log10 (energy [eV])", fontsize=22)
#   ax.set_ylabel(r"effective area [$km^{2}$]", fontsize=22)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   ax.set_xscale('log')
#   # ax.set_ylim(0,4)
#   ax.grid(True,alpha=0.2)
#   ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/effectiveArea.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # plotEffectiveArea(evtList,energyBins,weighting=False)
# # plotEffectiveArea(evtList_contained,energyBins,weighting=False)


# def plotFluxEnergyZenith(hdfFileList,energyBins,weighting):
#   '''
#   plots trigger efficiency in different zenith bins
#   '''
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   # colorIter = iter(colorsCustom)
#   colorIter = iter(colorsList)
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#     highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#     evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
#     energyList = []
#     sta3_fluxList = []
#     sta1_fluxList = []
#     totalEvts_list = []
#     for ebin, ebinStart in enumerate(energyBins[:-1]):
#       lowEdge_E = energyBins[ebin]
#       highEdge_E = energyBins[ebin+1]
#       evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
#       # totalEvts = len(evtEBin)
#       # weights = [ievt.H4aWeight for ievt in evtEBin]
#       weights = [ievt.directWeight for ievt in evtEBin]
#       totalEvts = np.sum(weights)
#       sta3_flux = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
#       sta1_flux = [ievt.STA1Trigger*ievt.H4aWeight for ievt in evtEBin]
#       # sta3_fluxList.append(sum(sta3_flux))
#       sta3_fluxList.append(sum(sta3_flux)/totalEvts*100.0)
#       totalEvts_list.append(totalEvts)
#       energyList.append((lowEdge_E+highEdge_E)/2.0)
#       # print("zen energy frac",binStart,np.log10(ebinStart),sum(sta3_flux)/totalEvts*100.0,sum(sta1_flux)/totalEvts*100.0)
#       # if sum(sta3_flux) > 10**-3:
#       #   ax.text((lowEdge_E+highEdge_E)/2.0,sum(sta3_flux),s="{0:.2f}".format(sum(sta3_flux)/totalEvts*100.0),fontsize=8)    
#     refLine=ax.plot(energyList,sta3_fluxList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
#     ax.plot(energyList,totalEvts_list,".",ls='-',lw = 2.5,c=refLine[0].get_color(),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"energy [GeV]", fontsize=22)
#   # ax.set_ylabel(r"rate [Hz]", fontsize=22))
#   ax.set_ylabel(r"percentage", fontsize=22)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   ax.set_xscale('log')
#   # ax.set_yscale('log')
#   ax.grid(True,alpha=0.2)
#   ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/fluxEnergyZenith.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # plotFluxEnergyZenith(hdf5NullList,energyBins,weighting=False)
# def plotEnergySpectraCompare(evtList,evtListOfficial,filtKey,filtKeyOfficial,yscale,suffix,energyScale):
#   if filtKey == "IceTopSTA5_13_filter":
#     evtListTrig = [ievt for ievt in evtList if abs(ievt.ITSTA5_filter-1)<0.01]
#   elif filtKey == "SDST_IceTopSTA3_13_filter":
#     evtListTrig = [ievt for ievt in evtList if abs(ievt.ITSTA3_filter-1)<0.01]
#   if filtKeyOfficial == "IceTopSTA5_12_filter":
#     evtListTrigOfficial = [ievt for ievt in evtListOfficial if abs(ievt.ITSTA5_filter-1)<0.01]
#   elif filtKeyOfficial == "SDST_IceTopSTA3_12_filter":
#     evtListTrigOfficial = [ievt for ievt in evtListOfficial if abs(ievt.ITSTA3_filter-1)<0.01]
#   else:
#     print("filtkey",filtKeyOfficial)
#   energy = np.log10([ievt.energy for ievt in evtListTrig])
#   weights = [ievt.H4aWeight for ievt in evtListTrig]
#   energyOfficial = np.log10([ievt.energy for ievt in evtListTrigOfficial])
#   weightsOfficial = [ievt.H4aWeight for ievt in evtListTrigOfficial]
#   # print("lengths",len(weights),weights,len(energy))
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   hitBins = np.linspace(5.0,8.0,31)
#   hist,binEdge = np.histogram(energy,hitBins,weights=weights)
#   histOfficial,binEdge = np.histogram(energyOfficial,hitBins,weights=weightsOfficial)
#   print("sum of official rate",sum(histOfficial))
#   print("sum of simulated rate",sum(hist))
#   ax.text(0.05, 0.95, r"{} official rate:{:.1f} Hz, sim. rate:{:.1f} Hz".format(filtKeyOfficial,sum(histOfficial),sum(hist)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
#   binCenter = [(binEdge[i+1]+binEdge[i])/2.0 for i,j in enumerate(binEdge[:-1])]
#   H = [h*(10**E)**energyScale for h,E in zip(hist,binCenter)]
#   HOfficial = [h*(10**E)**energyScale for h,E in zip(histOfficial,binCenter)]
#   ######################################################################################
#   # ax.hist(energy,bins=hitBins,histtype="step",lw=2.5,label=r"before weighting",alpha=1)
#   # ax.hist(energyOfficial,bins=hitBins,histtype="step",lw=2.5,label=r"after weighting",alpha=1)
#   ####################################################################################
#   ax.step(binCenter,H,"-",where="mid",lw=2.5,label="my simulation",alpha=0.6)
#   ax.step(binCenter,HOfficial,"-",where="mid",lw=2.5,label="official simulation",alpha=0.6)
#   ####################################################################################
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=weights,lw=2.5,label="after weighting",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"log10(energy [GeV])", fontsize=20)
#   ax.set_ylabel(r"$E^{0:.1f}$ rate [Hz]".format(energyScale), fontsize=20)
#   ax.set_yscale(yscale)
#   # ax.set_ylim(10**10,10**11.5)
#   # ax.set_ylim(10**6,10**12)
#   # ax.set_xlim(6.5,7.5)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.2)
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/energySpec"+str(suffix)+"Compare.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # plotEnergySpectraCompare(evtList,evtListOfficial,"IceTopSTA5_13_filter","IceTopSTA5_12_filter","log","fluxLogSTA5",energyScale=1.8)
# # plotEnergySpectraCompare(evtList,evtListOfficial,"IceTopSTA5_13_filter","IceTopSTA5_12_filter","log","fluxLogSTA5",energyScale=1.8)
# # plotEnergySpectraCompare(evtList,evtListOfficial,"SDST_IceTopSTA3_13_filter","SDST_IceTopSTA3_12_filter","log","fluxLogSTA3",energyScale=1.8)


# def delta_t_hist(eventList,histBins,suffix):
#   if "HLC" in suffix:
#     deltaTList = [ievt.deltaTHLC for ievt in eventList]
#   elif "SLC" in suffix:
#     deltaTList = [ievt.deltaT3SLC for ievt in eventList]
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   hitBins = np.linspace(min(deltaTList),max(deltaTList),200)
#   ax.hist(deltaTList,bins=histBins,histtype="step",lw=2.5,alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"$\Delta_t$ [ns]", fontsize=20)
#   ax.set_ylabel(r"count", fontsize=20)
#   ax.set_yscale('log')
#   ax.set_xlim(0,10000)
#   ax.grid(True,alpha=0.2)
#   ax.set_title(suffix,fontsize=12)
#   # ax.legend(fontsize=20)
#   plt.savefig(plotFolder+"/"+str(suffix)+"delta_t.pdf",transparent=False,bbox_inches='tight')
#   plt.close()
# # delta_t_hist(evtList,"HLCVEM")
# # delta_t_hist(evtList,np.linspace(0,10000,2000),"SLCVEM")



# def getEventsZenith_(hdfFile,zenLim):
#   '''get events belonging to given zenith bin:
#   zenith angles in degree
#   ''' 
#   df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
#   events_MC = df_MCPrimary["Event"].values
#   zeniths_MC = df_MCPrimary["zenith"].values
#   selEvents = []
#   for ievent,izenith in zip(events_MC,zeniths_MC):
#     if np.rad2deg(izenith) >= zenLim[0] and np.rad2deg(izenith) < zenLim[1]:
#       selEvents.append(ievent)
#   return list(set(selEvents))

# def getTriggeredEvents_(hdfFile,zenLim,energyLim,weighting):
#   """
#   returns total events sta1 and sta3 trigered events
#   """
#   STA1Trigger_df = pd.read_hdf(hdfFile,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_isSTA1")  
#   STA3Trigger_df = pd.read_hdf(hdfFile,key="ITSMTTriggered")
#   # selEvents = getEventsZenith_(hdfFile,zenLim)
#   selEvents = getEventsZenithEnergy_(hdfFile,zenLim,energyLim)
#   eventList = STA1Trigger_df["Event"].values
#   totalEvts = np.ones(len(eventList))
#   if weighting == True:
#     # weights,_ = weightCalc(hdfFile)
#     weights = getValue_(hdfFile,"H4aWeight")
#     totalEvts = np.multiply(totalEvts,weights)
#     STA1Triggers = np.multiply(STA1Trigger_df["value"].values,weights)
#     STA3Triggers = np.multiply(STA3Trigger_df["value"].values,weights)
#   else:
#     STA1Triggers = STA1Trigger_df["value"].values
#     STA3Triggers = STA3Trigger_df["value"].values
#   STA1Triggers_select = []
#   STA3Triggers_select = []
#   totalEvts_select = []
#   for ievt,ista1,ista3,isEvt in zip(eventList,STA1Triggers,STA3Triggers,totalEvts):
#     if ievt in selEvents:
#       STA1Triggers_select.append(ista1)
#       STA3Triggers_select.append(ista3)
#       totalEvts_select.append(isEvt)
#   return np.sum(totalEvts),np.sum(STA1Triggers_select),np.sum(STA3Triggers_select)
  
# def getTriggeredEvents(hdfFileList,zenLim,energyLim,weighting):
#   """
#   sums total sta1 and sta3 events from all hdfFileList files
#   """
#   totalEvts = 0
#   sta1_evts = 0
#   sta3_evts = 0
#   for ihdf in hdfFileList:
#     n_total,n_sta1,n_sta3 = getTriggeredEvents_(ihdf,zenLim,energyLim,weighting)
#     totalEvts += n_total
#     sta1_evts += n_sta1
#     sta3_evts += n_sta3
#   return totalEvts,sta1_evts,sta3_evts

# def getEventsZenithEnergy_(hdfFile,zenLim,energyLim):
#   '''get events belonging to given zenith bin:
#   zenith angles in degree
#   ''' 
#   df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
#   events_MC = df_MCPrimary["Event"].values
#   zeniths_MC = df_MCPrimary["zenith"].values
#   energy_MC = df_MCPrimary["energy"].values
#   selEvents = []
#   for ievent,izenith,ienergy in zip(events_MC,zeniths_MC,energy_MC):
#     if np.rad2deg(izenith) >= zenLim[0] and np.rad2deg(izenith) < zenLim[1]:
#       if ienergy >= energyLim[0] and ienergy < energyLim[1]:
#         selEvents.append(ievent)
#   return list(set(selEvents))


# def trigZen(hdfFileList,zenBins,fraction,weighting,suffix):
#   """
#   plots number of triggered events in different zenith angle bins
#   """
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   angleBins = []
#   sta1_list = []
#   sta3_list = []
#   totalEvts_list = []
#   sta1Frac_list = []
#   sta3Frac_list = []
#   for n,bin_start in enumerate(zenBins[:-1]):
#     lowEdge = zenBins[n]
#     highEdge = zenBins[n+1]
#     totalEvts,sta1_evts,sta3_evts = getTriggeredEvents(hdfFileList,[lowEdge,highEdge],[5.0,8.0],weighting)
#     angleBins.append((lowEdge+highEdge)/2.0)
#     sta1_list.append(sta1_evts)
#     sta3_list.append(sta3_evts)
#     totalEvts_list.append(totalEvts)
#   if fraction == True:
#     sta1Frac_list = np.asarray(sta1_list)/np.asarray(totalEvts_list)
#     sta3Frac_list = np.asarray(sta3_list)/np.asarray(totalEvts_list)
#     ax.plot(angleBins,sta1Frac_list,"o-",c=next(colorsIter),label="STA1",alpha=1)
#     ax.plot(angleBins,sta3Frac_list,"o-",c=next(colorsIter),label="STA3",alpha=1)
#     ax.set_ylabel(r"fraction", fontsize=20)
#   elif fraction == False:
#       ax.plot(angleBins,sta1_list,"o-",c=next(colorsIter),label="STA1",alpha=1)
#       ax.plot(angleBins,sta3_list,"o-",c=next(colorsIter),label="STA3",alpha=1)
#       # ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
#       ax.set_ylabel(r"rate [Hz]", fontsize=20)
#       for iangle,ista1,ista3,itot in zip(angleBins,sta1_list,sta3_list,totalEvts_list):
#         # ax.text(iangle,ista1+1,s="{0:.3f}".format(ista1/itot))
#         # ax.text(iangle,ista3+1,s="{0:.3f}".format(ista3/itot))
#         # ax.text(iangle,itot-200,s="{0:.0f}".format(itot))
#         ax.text(iangle,ista1,s="{0:.4f}%".format(ista1/itot*100))
#         ax.text(iangle,ista3,s="{0:.4f}%".format(ista3/itot*100))
#         # ax.text(iangle,itot,s="{0:.0f}".format(itot))
#   ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
#   ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
#   # ax.set_yscale('log')
#   ax.set_title("total rate: STA1 = {:.2f} Hz, STA3 = {:.2f} Hz".format(np.sum(sta1_list),np.sum(sta3_list)),fontsize=12)
#   ax.set_xlim(0,65)
#   ax.grid(True,alpha=0.2)
#   ax.legend(fontsize=14)
#   # ax.legend(fontsize=14,ncol=2)
#   plt.savefig(plotFolder+"/"+str(suffix)+"Frac"+str(fraction)+"Wts"+str(weighting)+"trigger_zen.pdf",transparent=False,bbox_inches='tight')
#   plt.close()


# # trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,False,"HLCVEM")
# # trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],False,False,"HLCVEM")
# # trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,True,"HLCVEM")
# # trigZen(hdf5NullList,[0,10,20,30,40,50,65.1],False,True,"HLCVEM")

# def plotFluxEnergyZenithScaled(hdfFileList,energyBins,weighting):
#   '''
#   plots trigger efficiency in different zenith bins
#   '''
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   # colorIter = iter(colorsCustom)
#   colorIter = iter(colorsList)
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#     highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#     evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
#     energyList = []
#     totalEvts_list = []
#     sta3_fluxList = []
#     sta1_fluxList = []
#     for ebin, ebinStart in enumerate(energyBins[:-1]):
#       lowEdge_E = energyBins[ebin]
#       highEdge_E = energyBins[ebin+1]
#       evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
#       # totalEvts = len(evtEBin)
#       weights = [ievt.H4aWeight for ievt in evtEBin]
#       totalEvts = np.sum(weights)
#       sta3_flux = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
#       sta3_fluxList.append(sum(sta3_flux)*((lowEdge_E+highEdge_E)/2.0)**1.8)
#       totalEvts_list.append(totalEvts*((lowEdge_E+highEdge_E)/2.0)**1.8)
#       energyList.append((lowEdge_E+highEdge_E)/2.0)   
#     refLine=ax.plot(energyList,sta3_fluxList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
#     ax.plot(energyList,totalEvts_list,".",ls='-',lw = 2.5,c=refLine[0].get_color(),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"energy [GeV]", fontsize=22)
#   ax.set_ylabel(r"$E^{1.8}$ rate [Hz]", fontsize=22)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   ax.set_xscale('log')
#   ax.set_yscale('log')
#   ax.grid(True,alpha=0.2)
#   # ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/fluxEnergyZenithScaled.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # plotFluxEnergyZenithScaled(hdf5NullList,energyBins,weighting=False)



# # def plotTrigEfficiency(hdfFileList,energyBins,weighting):
# #   '''
# #   plots trigger efficiency in different zenith bins
# #   '''
# #   fig = plt.figure(figsize=(8,5))
# #   gs = gridspec.GridSpec(nrows=1,ncols=1)
# #   ax = fig.add_subplot(gs[0])
# #   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
# #     lowEdge = sin2ZenBins[nbin]
# #     highEdge = sin2ZenBins[nbin+1]
# #     energyList = []
# #     sta3_efficiencyList = []
# #     sta1_efficiencyList = []
# #     for ebin, ebinStart in enumerate(energyBins[:-1]):
# #       lowEdge_E = energyBins[ebin]
# #       highEdge_E = energyBins[ebin+1]
# #       totalEvts,sta1_evts,sta3_evts = getTriggeredEvents(hdfFileList,[lowEdge,highEdge],[lowEdge_E,highEdge_E],weighting)
# #       sta3_trigEff = triggerEfficiency(sta3_evts,totalEvts)
# #       sta1_trigEff = triggerEfficiency(sta1_evts,totalEvts)
# #       sta3_efficiencyList.append(sta3_trigEff)
# #       sta1_efficiencyList.append(sta1_trigEff)
# #       energyList.append((lowEdge_E+highEdge_E)/2.0)   
# #     ax.plot(energyList,sta3_efficiencyList,".",ls='-',lw = 2.5,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
# #   ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
# #   ax.set_xlabel(r"energy [GeV]", fontsize=24)
# #   ax.set_ylabel(r"trigger efficiency", fontsize=24)
# #   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
# #   ax.set_xscale('log')
# #   ax.grid(True,alpha=0.2)
# #   ax.legend(fontsize=12)
# #   plt.savefig(plotFolder+"/trigEfficiency.pdf",transparent=False,bbox_inches='tight')
# #   plt.close()

# # plotTrigEfficiency(hdf5NullList,energyBins,weighting=False)

# def plotEffectiveArea(evtList,energyBins,weighting,triggerType,containment):
#   '''
#   plots trigger efficiency in different zenith bins
#   '''
#   if containment == True:
#     evtList = containedEvents(evtList,640)
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   # colorIter = iter(colorsCustom)
#   colorIter = iter(colorsList)
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
#     highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
#     evtZenBin = [ievt for ievt in evtList if lowEdge <= ievt.zenith < highEdge]
#     energyList = []
#     efficiencyList = []
#     for ebin, ebinStart in enumerate(energyBins[:-1]):
#       lowEdge_E = energyBins[ebin]
#       highEdge_E = energyBins[ebin+1]
#       evtEBin = [ievt for ievt in evtZenBin if lowEdge_E <= ievt.energy < highEdge_E]
#       # totalEvts = len(evtEBin)
#       weights = [ievt.H4aWeight for ievt in evtEBin]
#       totalEvts = len(evtEBin)
#       # sta3 = [ievt.ITSMTTriggered*ievt.H4aWeight for ievt in evtEBin]
#       if triggerType == "sta3":
#         triggerList = [ievt.ITSMTTriggered for ievt in evtEBin]
#       elif triggerType == "sta1":
#         triggerList = [ievt.STA1Trigger for ievt in evtEBin]
#       elif triggerType == "slc3":
#         triggerList = [ievt.slc3Trig for ievt in evtEBin]
#       elif triggerType == "slc4":
#         triggerList = [ievt.slc4Trig for ievt in evtEBin]
#       elif triggerType == "slc5":
#         triggerList = [ievt.slc5Trig for ievt in evtEBin]
#       effArea = effectiveArea(sum(triggerList),totalEvts,showerArea(np.log10(lowEdge_E)))
#       efficiencyList.append(effArea)
#       energyList.append(np.log10((lowEdge_E+highEdge_E)/2.0*10**9))   
#     ax.plot(energyList,efficiencyList,".",ls='-',lw = 2.5,c=next(colorIter),label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#   ax.set_xlabel(r"log10 (E [eV])", fontsize=22)
#   ax.set_ylabel(r"effective area [$km^{2}$]", fontsize=22)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   # ax.set_xscale('log')
#   ax.set_ylim(0,5)
#   ax.set_xlim(14,17)
#   # ax.yaxis.set_minor_locator(MultipleLocator(100))
#   ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#   ax.grid(True,alpha=0.5)
#   ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/trig"+str(triggerType)+"cont"+str(containment)+"EffectiveArea.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="sta1",containment=False)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="sta3",containment=False)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc3",containment=False)
# # # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc4",containment=False)
# # # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc5",containment=False)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="sta1",containment=True)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="sta3",containment=True)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc3",containment=True)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc4",containment=True)
# # plotEffectiveArea(evtList,energyBins,weighting=False,triggerType="slc5",containment=True)





# def plotEffectiveArea(zenChargeList):
#   '''
#   plots trigger efficiency in different zenith bins
#   '''
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(nrows=1,ncols=1)
#   ax = fig.add_subplot(gs[0])
#   for nbin, binStart in enumerate(sin2ZenBins[:-1]):
#     zclist = [zc for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
#     smtTrigEfficiencyList, globalTrigEfficiencyList = energyEfficiency(energyBins,zclist,"effArea")
#     ax.plot(energyBins[:-1],smtTrigEfficiencyList,".",ls='-',lw = 2.5,label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi),alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
#   ax.set_xlabel(r"energy [GeV]", fontsize=24)
#   ax.set_ylabel(r"effective area [$km^{2}$]", fontsize=24)
#   # ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
#   ax.set_xscale('log')
#   ax.grid(True,alpha=0.2)
#   ax.legend(fontsize=12)
#   plt.savefig(plotFolder+"/trigEffectiveArea.pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# def rateZen(hdfFileList,zenBins,weighting,suffix):
#   """
#   plots number of triggered events in different zenith angle bins
#   """
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   angleBins = []
#   sta1_list = []
#   sta3_list = []
#   totalEvts_list = []
#   sta1Frac_list = []
#   sta3Frac_list = []
#   for n,bin_start in enumerate(zenBins[:-1]):
#     lowEdge = zenBins[n]
#     highEdge = zenBins[n+1]
#     totalEvts,sta1_evts,sta3_evts = getTriggeredEvents(hdfFileList,[lowEdge,highEdge],weighting)
#     angleBins.append((lowEdge+highEdge)/2.0)
#     sta1_list.append(sta1_evts)
#     sta3_list.append(sta3_evts)
#     totalEvts_list.append(totalEvts)
#   ax.plot(angleBins,np.asarray(sta1_list)/trigWindow,"o-",c=next(colorsIter),label="STA1",alpha=1)
#   ax.plot(angleBins,np.asarray(sta3_list)/trigWindow,"o-",c=next(colorsIter),label="STA3",alpha=1)
#   # ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
#   ax.set_ylabel(r"rate [Hz]", fontsize=20)
#   ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
#   ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
#   # ax.set_yscale('log')
#   ax.set_xlim(0,65)
#   ax.set_title("total rate STA1 {:.2f} Hz STA3 {:.2f} Hz".format(np.sum(sta1_list)/trigWindow,np.sum(sta3_list)/trigWindow),fontsize=12)
#   ax.grid(True,alpha=0.2)
#   ax.legend(fontsize=14)
#   # ax.legend(fontsize=14,ncol=2)
#   plt.savefig(plotFolder+"/"+str(suffix)+"Wts"+str(weighting)+"triggerRate_zen.pdf",transparent=False,bbox_inches='tight')
#   plt.close()
# # rateZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,"HLCVEM")
# # rateZen(hdf5NullList,[0,10,20,30,40,50,65.1],True,"HLCVEM")





# def delta_t_events_(hdfFile,events,key):
#   '''gets delta_t out of given hdffile
#   belonging to given events
#   '''
#   df_delta_t = pd.read_hdf(hdfFile,key=key)
#   delta_t_selected = []
#   events_t = df_delta_t["Event"].values
#   delta_t = df_delta_t["item"].values
#   for jevent,jdelta_t in zip(events_t,delta_t):
#     if jevent in events:
#       delta_t_selected.append(jdelta_t)
#   return np.asarray(delta_t_selected)
# def delta_t_zenith_(hdfFile,zenLim,key):
#   sel_evts = getEventsZenith_(hdfFile,zenLim)
#   return delta_t_events_(hdfFile,sel_evts,key)

# def delta_t_zenith(hdfFileList,zenLim,key):
#   delta_t_list = []
#   for ihdf in hdfFileList:
#     delta_t = delta_t_zenith_(ihdf,zenLim,key)
#     delta_t_list = np.concatenate((delta_t_list,delta_t))
#   return delta_t_list

# def delta_t_hist_zen_bins(hdfFileList,zenBins,suffix,key):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     delta_ts = getVectorItem(hdfFileList,key)
#     lastBin = min(3000,max(delta_ts))
#     delta_t_bins = np.linspace(min(delta_ts),lastBin,300)
#     for n,bin_start in enumerate(zenBins[:-1]):
#       lowEdge = zenBins[n]
#       highEdge = zenBins[n+1]
#       delta_t_this = delta_t_zenith(hdfFileList,zenLim=[lowEdge,highEdge],key=key)
#       print("plotting hist",lowEdge,highEdge,len(delta_t_this))
#       ax.hist(delta_t_this,bins=delta_t_bins,
#         histtype="step",label=r"$\theta$ = {0:.0f}$^{{\circ}}$-{1:.0f}$^{{\circ}}$, {2:d} evts".format(lowEdge,highEdge,len(delta_t_this)),color=next(colorsIter),lw=2.5,alpha=1)
#     ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
#     ax.set_xlabel(r"$\Delta_t$ [ns]", fontsize=20)
#     ax.set_ylabel(r"count", fontsize=20)
#     ax.set_yscale('log')
#     # ax.set_xscale('log')
#     # ax.set_xlim(None,2000)
#     ax.grid(True,alpha=0.2)
#     ax.set_title(key,fontsize=12)
#     ax.legend(fontsize=14)
#     # ax.legend(fontsize=14,ncol=2)
#     plt.savefig(plotFolder+"/"+str(suffix)+"delta_t_zenBins.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # delta_t_hist_zen_bins(hdf5NullList,[0,10,20,30,40,50,65.1],"HLCVEM","OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_delta_t")
# # delta_t_hist_zen_bins(hdf5NullList,[0,10,20,30,40,50,65.1],"SLCVEM","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_delta_t")

# def getZenith_(hdfFile):
#   df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
#   return np.rad2deg(df_MCPrimary["zenith"].values)
# def getZenith(hdfFileList):
#   zenithList = np.array([])
#   for ihdf in hdfFileList:
#     izenith = getZenith_(ihdf)
#     print("izenith")
#     if len(izenith)>0:
#       zenithList = np.concatenate((zenithList,izenith))
#   return zenithList

# def getZenithBug_(hdfFile):
#   df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
#   df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
#   return np.rad2deg(df_MCPrimary["zenith"].values)
# def getZenith(hdfFileList):
#   zenithList = np.array([])
#   for ihdf in hdfFileList:
#     izenith = getZenith_(ihdf)
#     print("izenith")
#     if len(izenith)>0:
#       zenithList = np.concatenate((zenithList,izenith))
#   return zenithList

# def histZen(hdfFileList,suffix):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     zenithList = getZenith(hdfFileList)
#     zenithBins = np.linspace(min(zenithList),max(zenithList),100)
#     ax.hist(zenithList,bins=zenithBins,histtype="step",lw=2.5,alpha=1)
#     ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
#     ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
#     ax.set_ylabel(r"count", fontsize=20)
#     ax.set_yscale('log')
#     # ax.set_xlim(None,2000)
#     ax.grid(True,alpha=0.2)
#     # ax.set_title(fontsize=12)
#     ax.legend(fontsize=14)
#     plt.savefig(plotFolder+"/"+str(suffix)+"zenithDist.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # histZen(hdf5NullList,"hist")

# def applyWeight(hitList,weights):
#   return np.multiply(hitList,weights)

# SLCHits = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalHit")
# # SLCHitsW = applyWeight(SLCHits,adjustedWeights)
# SLCDuration = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration")
# SLC_Q = getValue(hdf5NullList,"OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalCharge")
# # SLC_QW = applyWeight(SLC_Q,adjustedWeights)

# HLCHits = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalHit")
# # HLCHitsW = applyWeight(HLCHits,adjustedWeights)
# HLCDuration = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeHitTimeDuration")
# HLC_Q = getValue(hdf5NullList,"OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalCharge")
# # HLC_QW = applyWeight(HLC_Q,adjustedWeights)




# def getTimeCharge(key,hdf5List):
#     timeList = np.array([])
#     chargeList = np.array([])
#     for ihdf in hdf5List:
#         time,charge = getTimeCharge_(key,ihdf)
#         timeList = np.concatenate((timeList,time))
#         chargeList = np.concatenate((chargeList,charge))
#     return timeList, chargeList

# def getTimeCharge_(key,hdfFile):
#     '''takes HLC/SLC pulses key and hdf file, returns time and charge
#     '''
#     dataT = pd.read_hdf(hdfFile,key=key)
#     return dataT["time"].values,dataT["charge"].values


# def plotTimeHistogram(time,suffix,key):
#     time = np.log10(time)
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     hitBins = np.linspace(min(time),max(time),200)
#     ax.hist(time,bins=hitBins,histtype="step",lw=2.5,alpha=1)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#     ax.set_xlabel(r"log10(time [ns])", fontsize=20)
#     ax.set_ylabel(r"count", fontsize=20)
#     ax.set_yscale('log')
#     ax.grid(True,alpha=0.2)
#     ax.set_title(key,fontsize=16)
#     # ax.legend(fontsize=20)
#     plt.savefig(plotFolder+"/"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
#     plt.close()

# def plotTimeChargeScatter(time,charge,suffix,key):
#     time = np.log10(time)
#     charge = np.log10(charge)
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     ax.scatter(time,charge,s=10)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#     ax.set_xlabel(r"log10(time [ns])", fontsize=20)
#     ax.set_ylabel(r"log10(charge)", fontsize=20)
# #     ax.set_yscale('log')
#     ax.set_ylim(-4.5,4.5)
#     ax.set_title(key,fontsize=16)
#     ax.grid(True,alpha=0.2)
#     plt.savefig(plotFolder+"/scatterChargeTime"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
#     plt.close()
    


# def scatter_hist_(x, y, ax, ax_histx, ax_histy):
#     # no labels
#     ax_histx.tick_params(axis="x", labelbottom=False)
#     ax_histy.tick_params(axis="y", labelleft=False)

#     # the scatter plot:
#     ax.scatter(x, y)

#     # now determine nice limits by hand:
#     binwidth = 0.25
#     xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
#     lim = (int(xymax/binwidth) + 1) * binwidth

#     bins = np.arange(-lim, lim + binwidth, binwidth)
#     ax_histx.hist(x, bins=bins)
#     ax_histy.hist(y, bins=bins, orientation='horizontal')

# def scatter_hist(x,y,suffix):
#     fig = plt.figure(figsize=(8, 5))
#     gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7),
#                           left=0.1, right=0.9, bottom=0.1, top=0.9,
#                           wspace=0.05, hspace=0.05)
#     ax = fig.add_subplot(gs[1, 0])
#     ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
#     ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
#     # use the previously defined function
#     scatter_hist_(x, y, ax, ax_histx, ax_histy)
#     plt.savefig(plotFolder+"/scatterChargeTimeHist"+str(suffix)+".pdf",transparent=False,bbox_inches='tight')
# #     plt.show()
#     plt.close()

# def plotScatterHist(hdfFileList,key,suffix):
#     timeList,chargeList = getTimeCharge(key=key,hdf5List=hdfFileList)
#     plotTimeHistogram(timeList,suffix,key)
#     plotTimeChargeScatter(timeList,chargeList,suffix,key)
# #     scatter_hist(timeList,chargeList,suffix)
# #     print("checking length",len(timeList),len(chargeList),len(hdfFileList))


# # plotScatterHist(hdf5NullList,key="OfflineIceTopSLCVEMPulses",suffix="SLCVEMTime")
# # plotScatterHist(hdf5NullList,key="OfflineIceTopHLCVEMPulses",suffix="HLCVEMTime")
# # plotScatterHist(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge",suffix="SLCVEMTimeCleaned")
# # plotScatterHist(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge",suffix="HLCVEMTimeCleaned")
# # dataT = pd.read_hdf(hdf5NullList[0],key="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration")


# def getTankHitDuration(keyDuration,keyHit,hdf5List):
#     hitDurationList = np.array([])
#     hitList = np.array([])
#     for ihdf in hdf5List:
#         ihitDuration,ihit = getTankHitDuration_(keyDuration,keyHit,ihdf)
#         hitDurationList = np.concatenate((hitDurationList,ihitDuration))
#         hitList = np.concatenate((hitList,ihit))
#     return hitDurationList,hitList

# def getTankHitDuration_(keyDuration,keyHit,hdfFile):
#     '''takes HLC/SLC pulses key and hdf file, returns time and charge
#     '''
#     dataDuration = pd.read_hdf(hdfFile,key=keyDuration)
#     dataHit = pd.read_hdf(hdfFile,key=keyHit)
#     return dataDuration["value"].values,dataHit["value"].values


# def plotHitDuration(hit,HitDuration,prefix):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     ax.scatter(hit,HitDuration,s=10)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#     ax.set_xlabel(r"tank hit", fontsize=20)
#     ax.set_ylabel(r"hit duration [ns]", fontsize=20)
#     # ax.set_yscale('log')
#     ax.set_ylim(0.1,None)
#     ax.grid(True,alpha=0.2)
#     plt.savefig(plotFolder+"/scatterHitDuration"+str(prefix)+".pdf",transparent=False,bbox_inches='tight')
#     plt.close()




# # slcHitDuration,slcHit = getTankHitDuration(keyDuration="OfflineIceTopSLCTankPulsesHitTimeDuration",keyHit="OfflineIceTopSLCTankPulsesTotalHit",hdf5List=hdf5NullList)
# # hlcHitDuration,hlcHit = getTankHitDuration(keyDuration="OfflineIceTopHLCTankPulsesHitTimeDuration",keyHit="OfflineIceTopHLCTankPulsesTotalHit",hdf5List=hdf5NullList)
# # plotHitDuration(slcHit,slcHitDuration,prefix="SLC")
# # plotHitDuration(hlcHit,hlcHitDuration,prefix="HLC")
# # slcHitDuration,slcHit = getTankHitDuration(keyDuration="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeHitTimeDuration",keyHit="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalHit",hdf5List=hdf5NullList)
# # hlcHitDuration,hlcHit = getTankHitDuration(keyDuration="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeHitTimeDuration",keyHit="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalHit",hdf5List=hdf5NullList)
# # plotHitDuration(slcHit,slcHitDuration,prefix="SLCClean")
# # plotHitDuration(hlcHit,hlcHitDuration,prefix="HLCClean")

# print("slc hit w",SLCHits)
# # print("slc hit w",SLCHitsW)
# # plotHitDuration(SLCHitsW,SLCDuration,"SLCWeightHit")


# def hitsWZenith_(hdfFile,key):
#   df_MCPrimary = pd.read_hdf(hdfFile,key="MCPrimary")
#   zenithList = np.rad2deg(df_MCPrimary["zenith"].values)
#   hitList = pd.read_hdf(hdfFile,key=key)["value"].values
#   return hitList,zenithList

# def hitsWZenith(hdfFileList,key):
#   hitList = np.array([])
#   zenithList = np.array([])
#   for ihdf in hdfFileList:
#     hits,zeniths = hitsWZenith_(ihdf,key)
#     hitList = np.concatenate((hitList,hits))
#     zenithList = np.concatenate((zenithList,zeniths))
#   return hitList,zenithList
# def plotScatterHitsZenith(hdfFileList,key,prefix):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(ncols=1,nrows=1)
#     ax = fig.add_subplot(gs[0])
#     hitList,zenithList = hitsWZenith(hdfFileList,key)
#     ax.scatter(zenithList,hitList,s=10)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#     ax.set_xlabel(r"$\theta ^{{\circ}}$", fontsize=20)
#     ax.set_ylabel(r"" + str(prefix)+" hits [ns]", fontsize=20)
#     # ax.set_yscale('log')
#     ax.set_title(key,fontsize=16)
#     ax.set_ylim(0.1,None)
#     ax.grid(True,alpha=0.2)
#     plt.savefig(plotFolder+"/scatterHitzenith"+str(prefix).replace(" ","")+".pdf",transparent=False,bbox_inches='tight')
#     plt.close()

# # plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopHLCVEMPulsesTotalHit",prefix="HLC VEM")
# # plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopSLCVEMPulsesTotalHit",prefix="SLC VEM")

# # plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalHit",prefix="Clean HLC VEM")
# # plotScatterHitsZenith(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalHit",prefix="Clean SLC VEM")

# def eventsZenith(zenLim):
#   dataT = pd.read_hdf(hdfFile,key="MCPrimary")
#   zeniths = dataT["zenith"].values


# def compareTwoWeights(evtList):
#   frameWeight = [ievt.H4aWeight for ievt in evtList]
#   calcWeight = [ievt.directWeight for ievt in evtList]
#   print("length of two weights",len(frameWeight),len(calcWeight))
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   # hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
#   hitBins = np.linspace(min(min(calcWeight),min(calcWeight)),max(max(calcWeight),max(calcWeight)),100)
#   ax.hist(calcWeight,bins=hitBins,histtype="step",lw=2.5,label=r"Weight(calc)",alpha=1)
#   ax.hist(frameWeight,bins=hitBins,histtype="step",lw=2.5,label="Weight(frame)",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
#   ax.set_xlabel(r"weight", fontsize=20)
#   ax.set_ylabel(r"count", fontsize=20)
#   ax.set_yscale("log")
#   # ax.set_ylim(None,10**5)
#   # ax.set_xscale('log')
#   ax.grid(True,alpha=0.2)
#   # ax.set_title(key,fontsize=16)
#   ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/weightCompare.pdf",transparent=False,bbox_inches='tight')
#   plt.close()
# # compareTwoWeights(evtList)

# def stationRate(hdfFileList,key,prefix,hitBins,xlabel):
#   stationHits = getValue(hdfFileList,key)
#   weights = getValue(hdfFileList,"H4aWeight")
#   triggerSTA3 = getValue(hdfFileList,"ITSMTTriggered")
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   # hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
#   ax.hist(stationHits,bins=hitBins,histtype="step",weights=[w*t for w,t in zip(weights,triggerSTA3)],lw=2.5,label=r"stations",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=12)
#   ax.tick_params(axis='both',which='major',length=5)
#   ax.tick_params(axis='both',which='minor',length=3)
#   ax.set_xlabel(xlabel, fontsize=12)
#   ax.set_ylabel(r"rate [Hz]", fontsize=12)
#   ax.set_yscale("log")
#   ticks = hitBins
#   ax.set_xticks(ticks[::10])
#   ax.set_xticks(ticks[::2],True)
#   ax.set_xlim(0,None)
#   # ax.set_xscale('log')
#   # ax.grid(True,alpha=0.2)
#   # ax.set_title(key,fontsize=16)
#   # ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/Rate"+str(prefix)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # stationRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_hitStations",prefix="HLCTankStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (HLC Tank)")
# # stationRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_hitStations",prefix="SLCTankStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (SLC Tank)")
# # stationRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanCharge_hitTanks",prefix="HLCTankTank",hitBins=np.linspace(0,160,161),xlabel=r"$N_{tanks}$ per event (HLC Tank)")
# # stationRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanCharge_hitTanks",prefix="SLCTankTank",hitBins=np.linspace(0,80,81),xlabel=r"$N_{tanks}$ per event (SLC Tank)")

# # stationRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_hitStations",prefix="HLCVEMStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (HLC VEM)")
# # stationRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_hitStations",prefix="SLCVEMStation",hitBins=np.linspace(0,80,81),xlabel=r"$N_{station}$ per event (SLC VEM)")
# # stationRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge_hitTanks",prefix="HLCVEMTank",hitBins=np.linspace(0,160,161),xlabel=r"$N_{tanks}$ per event (HLC VEM)")
# # stationRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge_hitTanks",prefix="SLCVEMTank",hitBins=np.linspace(0,80,81),xlabel=r"$N_{tanks}$ per event (SLC VEM)")




# def chargeRate(hdfFileList,key,prefix,hitBins,xlabel):
#   Qtot = np.log10(getValue(hdfFileList,key))
#   weights = getValue(hdfFileList,"H4aWeight")
#   triggerSTA3 = getValue(hdfFileList,"ITSMTTriggered")
#   fig = plt.figure(figsize=(8,5))
#   gs = gridspec.GridSpec(ncols=1,nrows=1)
#   ax = fig.add_subplot(gs[0])
#   # hitBins = np.linspace(min(min(calcWeight),min(frameWeight)),max(max(calcWeight),max(frameWeight)),100)
#   ax.hist(Qtot,bins=hitBins,histtype="step",weights=[w*t for w,t in zip(weights,triggerSTA3)],lw=2.5,label=r"stations",alpha=1)
#   # ax.hist(energy,bins=hitBins,histtype="step",weights=adjWeights,lw=2.5,label="after adj weighting",alpha=1)
#   ax.tick_params(axis='both',which='both', direction='in', labelsize=12)
#   ax.tick_params(axis='both',which='major',length=5)
#   ax.tick_params(axis='both',which='minor',length=3)
#   ax.set_xlabel(xlabel, fontsize=12)
#   ax.set_ylabel(r"rate [Hz]", fontsize=12)
#   ax.set_yscale("log")
#   ticks = hitBins
#   ax.set_xticks(ticks[::10])
#   ax.set_xticks(ticks[::2],True)
#   # ax.set_xlim(0,None)
#   # ax.set_xscale('log')
#   # ax.grid(True,alpha=0.2)
#   # ax.set_title(key,fontsize=16)
#   # ax.legend(fontsize=10,ncol=3,loc="lower center")
#   plt.savefig(plotFolder+"/chargeRate"+str(prefix)+".pdf",transparent=False,bbox_inches='tight')
#   plt.close()

# # chargeRate(hdf5NullList,key="OfflineIceTopHLCTankPulsesCleanTimeCleanChargeTotalCharge",prefix="HLCTankCharge",hitBins=np.linspace(0,5,81),xlabel=r"log($Q_{tot}/VEM$) per event (HLC tank)")
# # chargeRate(hdf5NullList,key="OfflineIceTopSLCTankPulsesCleanTimeCleanChargeTotalCharge",prefix="SLCTankCharge",hitBins=np.linspace(-2,4,81),xlabel=r"log($Q_{tot}/VEM$) per event (SLC tank)")
# # chargeRate(hdf5NullList,key="OfflineIceTopHLCVEMPulsesCleanTimeCleanChargeTotalCharge",prefix="HLCVEMCharge",hitBins=np.linspace(0,5,81),xlabel=r"log($Q_{tot}/VEM$) per event (HLC VEM)")
# # chargeRate(hdf5NullList,key="OfflineIceTopSLCVEMPulsesCleanTimeCleanChargeTotalCharge",prefix="SLCVEMCharge",hitBins=np.linspace(-2,4,81),xlabel=r"log($Q_{tot}/VEM$) per event (SLC VEM)")





