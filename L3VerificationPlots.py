#!/usr/bin/env python3

from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
from icecube.recclasses import I3LaputopParams, LaputopParameter as Par

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import glob
import numpy as np
import pandas as pd
import pickle
import h5py


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input',"-i", type=str, nargs='+', 
 default=["/home/enpaudel/icecube/triggerStudy/simFiles/ITGenTest/FeDAT000001GenDet.i3.bz2"], help='Input files for verification plots')
parser.add_argument("--isSim",action="store_true",
  default=False, dest="isSim",help="tests if simulation/data")
parser.add_argument("--writeh5",action="store_true",
  default=False, dest="writeh5",help="write h5 files")
args = parser.parse_args()


exceptionTanks_HG = {39:62,26:62,67:64,74:62}
exceptionTanks_LG = {26:61,67:63}

#example l3 data /data/ana/CosmicRay/IceTop_level3/exp/IC86.2020/2021/0305/Run00135057/Level3_IC86.2020_data_Run00135057_Subrun00000000_00000270.i3.bz2

'''
ex1
python L3VerificationPlots.py -i 
/data/ana/CosmicRay/IceTop_level3/exp/IC86.2020/2021/0301/Level3_IC86.2020_data_Run00135045_Subrun00000000_00000000.i3.bz2
ex2
python L3VerificationPlots.py -i /data/sim/IceTop/2023/generated/untriggered/level3/pDAT005898L3.i3.gz --isSim

ex3 just for the plot
python L3VerificationPlots.py
'''


plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"
############################ if need to check GRL#
GRL = "/data/ana/CosmicRay/IceTop_GRL/"
GRL2020 = GRL + "/IC86_2020_GoodRunInfo_4IceTop_level2.txt"
GRL2020_df = pd.read_csv(GRL2020,sep="\t", header=0,escapechar='#')
GRL_IT_List = np.asarray(GRL2020_df["OutDir"][GRL2020_df["Good_it"]==1])
GRL_RunNum = np.asarray(GRL2020_df["RunNum"][GRL2020_df["Good_it"]==1])
GRL_LiveTime = np.asarray(GRL2020_df["LiveTime"][GRL2020_df["Good_it"]==1])
# print("GRL_RunNum",GRL_RunNum)

Level3_2021_path = "/data/ana/CosmicRay/IceTop_level3/exp/IC86.2020/2021/"
RunListFiles = sorted(glob.glob(Level3_2021_path+"03??/*/*.i3.*"))
RunList = sorted(glob.glob(Level3_2021_path+"03??/*"))
RunList = [int(irun.split("/")[-1].replace("Run00","")) for irun in RunList]
RunList = [irun for irun in RunList if irun in GRL_RunNum]
LiveTimeList = [itime for irun,itime in zip(GRL_RunNum,GRL_LiveTime) if irun in RunList]
totalTime = sum(LiveTimeList)

# print("runlist",len(RunList),RunList)
# print("timelist",len(LiveTimeList),sum(LiveTimeList)/3600.0,LiveTimeList)
###############################################
#for simulation
Level3_2021_sim_path = "/data/sim/IceTop/2023/generated/untriggered/level3"
SimListFiles = sorted(glob.glob(Level3_2021_sim_path+"/*.i3.gz"))


###############################################

RunListFiles = sorted(glob.glob(Level3_2021_path+"03??/*/*.i3.*"))
liveTime = glob.glob(Level3_2021_path+"03??/*/*.pickle")[0]
# with (open(liveTime, "rb")) as timefile:
#   fobj = pickle.load(timefile)
# print("fobj",fobj)

# RunList = sorted(glob.glob(Level3_2021_path+"03??/*"))
# print("runlist",len(RunList),RunList)

#############################################################
hdf5Path = "/home/enpaudel/icecube/triggerStudy/plots/hdf5/"
if args.isSim:
  inputFileName = hdf5Path+"sim_"+args.input[0].split("/")[-1].split(".")[0]+".hdf5"
else:
  inputFileName = hdf5Path+args.input[0].split("/")[-1].split(".")[1]+".hdf5"
print("input",inputFileName)

# print(SimListFiles)
hdf5FileSimu = glob.glob(hdf5Path+"sim_*.hdf5")
hdf5FileData = glob.glob(hdf5Path+"2020_data*.hdf5")



histDict = {
  "nTanksHLC":[np.linspace(-1,162,164),[0,163],[10**(-6.7),10**0.5],"# of tanks per event","OfflineIceTopHLCTankPulses"],
  "nTanksSLC":[np.linspace(-1,100,102),[0,90],[10**(-6.7),10**0.5],"# of tanks per event","OfflineIceTopSLCTankPulses"],
  "nTanksHLCRT":[np.linspace(-1,162,164),[0,163],[10**(-6.7),10**0.5],"# of tanks per event","IceTopHLCSeedRTPulses"],
  "nTanksHLCLaputop":[np.linspace(-1,162,164),[0,163],[10**(-6.7),10**0.5],"# of tanks per event","IceTopLaputopSeededSelectedHLC"],
  "nTanksSLCLaputop":[np.linspace(-1,100,102),[0,60],[10**(-6.7),10**0.5],"# of tanks per event","IceTopLaputopSeededSelectedSLC"],
  "nStationsHLC":[np.linspace(-1,100,102),[0,83],[10**(-6.7),10**0.5],"# of stations per event","OfflineIceTopHLCTankPulses"],
  "nStationsSLC":[np.linspace(-1,100,102),[0,83],[10**(-6.7),10**0.5],"# of tanks per event","OfflineIceTopSLCTankPulses"],
  "nStationsHLCRT":[np.linspace(-1,100,102),[0,83],[10**(-6.7),10**0.5],"# of stations per event","IceTopHLCSeedRTPulses"],
  "nStationsHLCLaputop":[np.linspace(-1,100,102),[0,83],[10**(-6.7),10**0.5],"# of stations per event","IceTopLaputopSeededSelectedHLC"],
  "nStationsSLCLaputop":[np.linspace(-1,100,102),[0,83],[10**(-6.7),10**0.5],"# of stations per event","IceTopLaputopSeededSelectedSLC"],
  "QtotHLC":[np.linspace(-2,5,101),[-1,4.5],[10**(-6.7),10**0.5],r"log$_{10}$(Q$_{tot}$/VEM)","OfflineIceTopHLCTankPulses"],
  "QtotHLCRT":[np.linspace(-2,5,101),[-1,4.5],[10**(-6.7),10**0.5],r"log$_{10}$(Q$_{tot}$/VEM)","IceTopHLCSeedRTPulses"],
  "QtotHLCLaputop":[np.linspace(-2,5,101),[-1,4.5],[10**(-6.7),10**0.5],r"log$_{10}$(Q$_{tot}$/VEM)","IceTopLaputopSeededSelectedHLC"],
  "QtotSLCLaputop":[np.linspace(-2,5,101),[-1.5,3.5],[10**(-6.7),10**0.5],r"log$_{10}$(Q$_{tot}$/VEM)","IceTopLaputopSeededSelectedSLC"],
  "QtotSLC":[np.linspace(-2,5,101),[-1.5,3.5],[10**(-6.7),10**0.5],r"log$_{10}$(Q$_{tot}$/VEM)","OfflineIceTopSLCTankPulses"],
  "logS125":[np.linspace(-4,4,101),[-3.2,4],[10**(-6.7),10**0.5],r"log$_{10}$(S$_{125}$/VEM)",""],
  "beta":[np.linspace(-2,12,101),[-0.5,11],[10**(-6.7),10**0.5],r"$\beta$",""],
  }


class VerificationPlots(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
    self.AddParameter("simulation","True if simulation, False if Data",False)

  def Configure(self):
    self.simulation = self.GetParameter("simulation")
    self.fig = plt.figure(figsize=(8,5))
    self.gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(self.gs[0])
    self.nTanksHLC = []
    self.nTanksSLC = []
    self.nTanksHLCRT = []
    self.nTanksHLCLaputop = []
    self.nTanksSLCLaputop = []
    self.nStationsHLC = []
    self.nStationsSLC = []
    self.nStationsHLCRT = []
    self.nStationsHLCLaputop = []
    self.nStationsSLCLaputop = []
    self.QtotHLC = []
    self.QtotHLCRT = []
    self.QtotHLCLaputop = []
    self.QtotSLCLaputop = []
    self.QtotSLC = []
    self.logS125 = []
    self.beta = []


  def get_nstation(self,frame,pulseseries):
    '''
    returns number of tank hits and station hits

    '''
    hit_stations = []
    hit_doms = []
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
    for idom,pulses in psm:
      for pulse in pulses:
        hit_stations.append(idom.string)
        if ( (idom.om == (61 + (idom.string==26 or idom.string==39 or idom.string==74))) or (idom.om == (63 + (idom.string==67))) ):
          hit_doms.append(idom)
        break
    hit_stations = list(set(hit_stations))
    hit_tanks = list(set(hit_doms))
    return len(hit_tanks),len(hit_stations)

  def getQtot(self,frame,pulseseries):
    '''
    total charge in an event
    '''
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
    return sum([sum([ps.charge for ps in psm[i]]) for i in psm.keys()])


  def Physics(self,frame):
    """what if there is only one trigger in trigger hierarchy
    """
    # trigger = frame["I3Triggers"]
    # print(frame.keys())
    if frame.Has("H4aWeight"):
      iweight = frame["H4aWeight"].value
    else:
      iweight = 1.0
    if frame["I3EventHeader"].sub_event_stream == "IceTopSplit":
      if frame.Has("OfflineIceTopHLCTankPulses"):        
        nTanksHLC, nStationsHLC = self.get_nstation(frame,"OfflineIceTopHLCTankPulses")
        self.nTanksHLC.append([nTanksHLC,iweight])
        self.nStationsHLC.append([nStationsHLC,iweight])
        self.QtotHLC.append([np.log10(self.getQtot(frame,"OfflineIceTopHLCTankPulses")),iweight])
      if frame.Has("IceTopHLCSeedRTPulses"):
        nTanksHLCRT, nStationsHLCRT = self.get_nstation(frame,"IceTopHLCSeedRTPulses")
        self.nTanksHLCRT.append([nTanksHLCRT,iweight])
        self.nStationsHLCRT.append([nStationsHLCRT,iweight])
        self.QtotHLCRT.append([np.log10(self.getQtot(frame,"IceTopHLCSeedRTPulses")),iweight])
      if frame.Has("IceTopLaputopSeededSelectedHLC"):
        nTanksHLCLaputop, nStationsHLCLaputop = self.get_nstation(frame,"IceTopLaputopSeededSelectedHLC")
        self.nTanksHLCLaputop.append([nTanksHLCLaputop,iweight])
        self.nStationsHLCLaputop.append([nStationsHLCLaputop,iweight])
        self.QtotHLCLaputop.append([np.log10(self.getQtot(frame,"IceTopLaputopSeededSelectedHLC")),iweight])
      if frame.Has("IceTopLaputopSeededSelectedSLC"):
        nTanksSLCLaputop, nStationsSLCLaputop = self.get_nstation(frame,"IceTopLaputopSeededSelectedSLC")
        self.nTanksSLCLaputop.append([nTanksSLCLaputop,iweight])
        self.nStationsSLCLaputop.append([nStationsSLCLaputop,iweight])
        self.QtotSLCLaputop.append([np.log10(self.getQtot(frame,"IceTopLaputopSeededSelectedSLC")),iweight])
      if frame.Has("OfflineIceTopSLCTankPulses"):
        nTanksSLC, nStationsSLC = self.get_nstation(frame,"OfflineIceTopSLCTankPulses")
        self.nStationsSLC.append([nStationsSLC,iweight])
        self.nTanksSLC.append([nTanksSLC,iweight])
        self.QtotSLC.append([np.log10(self.getQtot(frame,"OfflineIceTopSLCTankPulses")),iweight])
      # if frame.Has("LaputopParams") and frame["Laputop"].fit_status == "OK":
      if frame.Has("LaputopParams") and str(frame["Laputop"].fit_status) == "OK":
        # print(type(frame["Laputop"].fit_status),str(frame["Laputop"].fit_status),str(frame["Laputop"].fit_status) == "OK")
        params = I3LaputopParams.from_frame(frame,"LaputopParams")
        self.logS125.append([params.value(Par.Log10_S125),iweight])
        self.beta.append([params.value(Par.Beta),iweight])
    self.PushFrame(frame)

  def plotHist(self,x,weights,xtype,bins):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(-1,max(nHits),max(nHits)+2)
    counts, bin_edges = np.histogram(x, bins=bins)
    # ax.hist(bin_edges[:-1],bins=bin_edges,weights=counts*1/totalTime,histtype="step",label="total",lw=1.5)
    ax.hist(x,bins=bins,weights=weights,histtype="step",label="total",lw=1.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel("{} per event".format(xtype), fontsize=22)
    ax.set_ylabel("rate [Hz]", fontsize=22)
    ax.set_yscale("log")
    ax.grid(True,alpha=0.4)
    # ax.text(0.82,0.88,s="mean {0} hits\n{1:.1f} per event".format(pulseDict[hitType],meanHit),size=15,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    # ax.legend(loc="upper left",fontsize=12)
    # ax.legend(fontsize=12)
    # ax.hist(weightspy3,histtype="step")
    plt.savefig(plotFolder+"/L3verification"+str(xtype)+".png",transparent=False,bbox_inches='tight')
    plt.close()

  def Finish(self):
    with h5py.File('{}'.format(inputFileName), 'w') as hf:
      for ix in histDict:
        print("key",ix)
        attr = getattr(self,ix)
        hf.create_dataset("{}".format(ix), data=attr)
        # values = [iattr[0] for iattr in attr]
        # weights = [iattr[1] for iattr in attr]
        # print("attr",attr,min(attr),max(attr))
        # self.plotHist(values,weights,ix,histDict[ix][0])
    # print("stations",self.nTanksHLC,self.nStationsHLC)
    # print("number of SMT triggered events",self.SMTTriggered )
    # print("number of SMT untriggered events",self.SMTUntriggered )
    # print("number of Global triggered events",self.GlobalTriggered)
    # print("number of Gloabl untriggered events",self.GlobalUntriggered)
    # hdftable = hdfwriter.I3HDFTableService(str(outputDir)+str(fileName)+"TrigCount.hdf5")

def plotHist(xdata,weightsdata,xsim,weightssim,xtype,bins,xlim,ylim,xlabel,text_label):
  fig = plt.figure(figsize=(8,5))
  gs = gridspec.GridSpec(nrows=1,ncols=1)
  ax = fig.add_subplot(gs[0])
  # bins = np.linspace(-1,max(nHits),max(nHits)+2)
  # counts, bin_edges = np.histogram(x, bins=bins)
  # ax.hist(bin_edges[:-1],bins=bin_edges,weights=counts*1/totalTime,histtype="step",label="total",lw=1.5)
  ax.hist(xdata,bins=bins,weights=weightsdata,histtype="step",label="data",lw=2)
  ax.hist(xsim,bins=bins,weights=weightssim,histtype="step",label="sim",lw=2)
  ax.tick_params(axis='both',which='both', direction='in', labelsize=18)
  # ax.set_xlabel("{} per event".format(xtype), fontsize=22)
  ax.set_xlabel("{}".format(xlabel), fontsize=18)
  ax.set_ylabel("rate [Hz]", fontsize=18)
  ax.text(.1, .98, text_label, ha='left', va='top', transform=ax.transAxes,fontsize=16)
  ax.set_xlim(xlim[0],xlim[1])
  ax.set_ylim(ylim[0],ylim[1])
  ax.set_yscale("log")
  ax.grid(True,alpha=0.4)
  # ax.text(0.82,0.88,s="mean {0} hits\n{1:.1f} per event".format(pulseDict[hitType],meanHit),size=15,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
  # ax.legend(loc="upper left",fontsize=12)
  ax.legend(fontsize=16)
  # ax.hist(weightspy3,histtype="step")
  plt.savefig(plotFolder+"/L3verification"+str(xtype)+".png",transparent=False,bbox_inches='tight')
  plt.close()

def getValueWeight(xvalue,hdf5List,isData):
  values = []
  weights = []
  for ihdf in hdf5List:
    with h5py.File(ihdf, 'r') as hf:
      arr = np.asarray(hf[xvalue][:]).T
      values += list(arr[0])
      if isData:
        weights += [dt/totalTime for dt in arr[1]]
      else:
        weights += list(arr[1])
  return values, weights




if args.writeh5 == True:
  tray = icetray.I3Tray()
  tray.AddModule("I3Reader","reader",
               # filenameList=RunListFiles,
               # filenameList=SimListFiles,
               filenameList=args.input,
              )

  tray.AddModule(VerificationPlots,"veriP",
              simulation=args.isSim,
              )
  # tray.AddModule("I3Writer","i3writer",
  #             filename=args.output,
  #             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
  #             )

  tray.Execute()
  tray.Finish()
else:
  for xvalue in histDict:
    print("plotting {}".format(xvalue))
    values_data,weights_data = getValueWeight(xvalue,hdf5FileData,isData=True)
    values_sim,weights_sim = getValueWeight(xvalue,hdf5FileSimu,isData=False)
    plotHist(values_data,weights_data,values_sim,weights_sim,xvalue,histDict[xvalue][0],
      histDict[xvalue][1],histDict[xvalue][2],histDict[xvalue][3],histDict[xvalue][4])




