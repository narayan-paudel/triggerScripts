#!/usr/bin/env python

import tables
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap

import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetFiltProcTrigCount.hdf5", help='Input files after running plotTrigHist.py.')
args = parser.parse_args()

f = tables.open_file("../simFiles/trigChkTest/DAT000005GenDetFiltProcTrigCount.hdf5")
print(f.root.OfflineIceTopSLCTankPulsesTotalTankHit)
# print(f.root.)

# f = tables.open_file("/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetFiltProcTrigCount.hdf5")
# print(f.root)
# print(f.root.ITTotalChargeSLC)
# print(f.root.ITTotalChargeSLC.cols.value[:])

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

# fh_SLCTank = pd.read_hdf("/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetFiltProcTrigCount.hdf5",key="OfflineIceTopSLCTankPulsesTotalCharge")
# fh_HLCTank = pd.read_hdf("/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetFiltProcTrigCount.hdf5",key="OfflineIceTopHLCTankPulsesTotalCharge")
fh_primary = pd.read_hdf(str(args.input[0]),key="MCPrimary")
zenithList = fh_primary['zenith'].values
# print("primary",fh_primary.columns)
print("primary",fh_primary['zenith'])
fh_SLCTank = pd.read_hdf(str(args.input[0]),key="OfflineIceTopSLCTankPulsesTotalCharge")
fh_HLCTank = pd.read_hdf(str(args.input[0]),key="OfflineIceTopHLCTankPulsesTotalCharge")

fh_SLCTankHit = pd.read_hdf(str(args.input[0]),key="OfflineIceTopSLCTankPulsesTotalTankHit")
fh_HLCTankHit = pd.read_hdf(str(args.input[0]),key="OfflineIceTopHLCTankPulsesTotalTankHit")

print(fh_HLCTankHit.columns)
# print(fh.head())
# print(fh.columns)
# print(fh["value"].values)
ITTotalChargeHLC = fh_HLCTank["value"].values
ITTotalChargeSLC = fh_SLCTank["value"].values

ITTotalHitHLC = fh_HLCTankHit["value"].values
ITTotalHitSLC = fh_SLCTankHit["value"].values
print("length compare",len(zenithList),len(ITTotalChargeHLC))

sin2ZenBins = [0.0,0.2,0.4,0.6,0.8,0.1]
class zenChargeHLCSLC(object):
	"""docstring for zenChargeHLCSLC"""
	def __init__(self, zenith,chargeHLC, chargeSLC,tankHitHLC,tankHitSLC):
		super(zenChargeHLCSLC, self).__init__()
		self.zenith = zenith
		self.chargeHLC = chargeHLC
		self.chargeSLC = chargeSLC
		self.tankHitHLC = tankHitHLC
		self.tankHitSLC = tankHitSLC

zenChargeList = []
for zen,hlc,slc,hlcHit,slcHit in zip(zenithList,ITTotalChargeHLC,ITTotalChargeSLC,ITTotalHitHLC,ITTotalHitSLC):
	zc = zenChargeHLCSLC(zen,hlc,slc,hlcHit,slcHit)
	zenChargeList.append(zc)

slcLim = 100
zenChargeList = np.asarray([zc for zc in zenChargeList if zc.chargeSLC < slcLim]) #to remove events with very high slc charge

def plotHistChargeZenBin(zenChargeList,LCType):
	'''
	plots histogram of charge for zenith bins.
	'''
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	if LCType == "HLC":
		hlcChargeList = sorted([zc.chargeHLC for zc in zenChargeList])
		chargeBins = np.linspace(min(hlcChargeList),max(hlcChargeList),40)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			chargeInBin = [zc.chargeHLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			ax.hist(chargeInBin ,bins=chargeBins,histtype="step",lw=2.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
	elif LCType == "SLC":
		slcChargeList = sorted([zc.chargeSLC for zc in zenChargeList])
		chargeBins = np.linspace(min(slcChargeList),max(slcChargeList),40)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			chargeInBin = [zc.chargeSLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			ax.hist(chargeInBin ,bins=chargeBins,histtype="step",lw=2.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"q [vem]", fontsize=24)
	ax.set_ylabel(r"$count$", fontsize=24)
	ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	ax.set_yscale('log')
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histChargeZenBins"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotHistChargeZenBin(zenChargeList,"HLC")
# plotHistChargeZenBin(zenChargeList,"SLC")


def plot2DHistChargeZenBin(zenChargeList,LCType):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	# zenithBins = [0.0,0.2,0.4,0.6,0.8,1.0]
	zenithBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
	zenithList = sorted([np.square(np.sin(zc.zenith)) for zc in zenChargeList])
	# zenithBins = np.linspace(min(zenithList),max(zenithList),4)
	if LCType == "HLC":
		hlcChargeList = sorted([zc.chargeHLC for zc in zenChargeList])
		hlcChargeBins = np.linspace(min(hlcChargeList),max(hlcChargeList),40)
		hist = ax.hist2d(zenithList,hlcChargeList,bins=(zenithBins,hlcChargeBins),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)
	elif LCType == "SLC":
		slcChargeList = sorted([zc.chargeSLC for zc in zenChargeList])
		slcChargeBins = np.linspace(min(slcChargeList),max(slcChargeList),40)
		hist = ax.hist2d(zenithList,slcChargeList,bins=(zenithBins,slcChargeBins),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)
	fig.colorbar(hist[3],ax=ax)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_ylabel(r"q [vem]", fontsize=24)
	ax.set_xlabel(r"$\sin ^2 \theta$", fontsize=24)
	ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/2DhistChargeZenBins"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plot2DHistChargeZenBin(zenChargeList,"HLC")
# plot2DHistChargeZenBin(zenChargeList,"SLC")
		
###################################################################################################################
#################################3for hit counts###################################################################
def plotHistHitZenBin(zenChargeList,LCType):
	'''
	plots histogram of charge for zenith bins.
	'''
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	if LCType == "HLC":
		hlcHitList = sorted([zc.tankHitHLC for zc in zenChargeList])
		hitBins = np.linspace(min(hlcHitList),max(hlcHitList),40)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			hitInBin = [zc.tankHitHLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			ax.hist(hitInBin ,bins=hitBins,histtype="step",lw=2.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
	elif LCType == "SLC":
		slcHitList = sorted([zc.tankHitSLC for zc in zenChargeList])
		hitBins = np.linspace(min(slcHitList),max(slcHitList),40)
		for nbin,binStart in enumerate(sin2ZenBins[:-1]):
			hitInBin = [zc.tankHitSLC for zc in zenChargeList if (np.square(np.sin(zc.zenith)) >= sin2ZenBins[nbin] and np.square(np.sin(zc.zenith)) < sin2ZenBins[nbin+1])]
			ax.hist(hitInBin ,bins=hitBins,histtype="step",lw=2.5,label=r"$\sin ^2 \theta$ = "+str(sin2ZenBins[nbin])+"-"+str(sin2ZenBins[nbin+1]),alpha=1)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_xlabel(r"q [vem]", fontsize=24)
	ax.set_ylabel(r"$count$", fontsize=24)
	ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	ax.set_yscale('log')
	ax.grid(True,alpha=0.2)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histHitsZenBins"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotHistHitZenBin(zenChargeList,"HLC")
plotHistHitZenBin(zenChargeList,"SLC")


def plot2DHistHitZenBin(zenChargeList,LCType):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	# zenithBins = [0.0,0.2,0.4,0.6,0.8,1.0]
	zenithBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
	zenithList = sorted([np.square(np.sin(zc.zenith)) for zc in zenChargeList])
	# zenithBins = np.linspace(min(zenithList),max(zenithList),4)
	if LCType == "HLC":
		hlcHitList = sorted([zc.tankHitHLC for zc in zenChargeList])
		hlcHitBins = np.linspace(min(hlcHitList),max(hlcHitList),40)
		hist = ax.hist2d(zenithList,hlcHitList,bins=(zenithBins,hlcHitBins),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)
	elif LCType == "SLC":
		slcHitList = sorted([zc.tankHitSLC for zc in zenChargeList])
		slcHitBins = np.linspace(min(slcHitList),max(slcHitList),40)
		hist = ax.hist2d(zenithList,slcHitList,bins=(zenithBins,slcHitBins),norm=mpl.colors.LogNorm(),cmin=None,cmax=None,alpha=1)
	fig.colorbar(hist[3],ax=ax)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=24)
	ax.set_ylabel(r"q [vem]", fontsize=24)
	ax.set_xlabel(r"$\sin ^2 \theta$", fontsize=24)
	ax.set_title("OfflineIceTop"+LCType+"TankPulses",fontsize=24)
	# ax.set_xscale('log')
	ax.grid(True,alpha=0.2)
	# ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/2DhistHitZenBins"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

plot2DHistHitZenBin(zenChargeList,"HLC")
plot2DHistHitZenBin(zenChargeList,"SLC")







def plotHistCharge(charges,LCType):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	chargeMean = np.mean(charges)
	chargeSD = np.std(charges)
	minLim = min(charges)
	maxLim = max(charges)
	n_bins = np.linspace(minLim,maxLim,200)
	# n_bins = np.linspace(0,90,9)
	ax.hist(charges,bins=n_bins,histtype="step",lw=2.5,edgecolor='#2ca02c',label=str(LCType)+" hist",alpha=1)
	ax.axvline(x=chargeMean, ymin=0, ymax=1,c='#2ca02c',ls="-.",label=r"$\bar{q}$")
	ax.axvspan(chargeMean - chargeSD,chargeMean + chargeSD,ymin=0, ymax=1,color='#2ca02c',ls="-",label=r"$\sigma$",alpha=0.1)
	# ax.hist(thetaList,bins=n_bins,histtype="step",lw=2.5,fc='orange',label="ratio: 12 - 75",alpha=0.75)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
	ax.tick_params(axis='both',which='major', length=8,width=1.2)
	ax.tick_params(axis='both',which='minor', length=4)
	ax.set_yscale('log')
	ax.grid(True)
	# ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
	# ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
	ax.set_xlabel(r"q [vem]", fontsize=28)
	ax.set_ylabel(r"$count$", fontsize=22)
	ax.set_ylim(1,None)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histCharge"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()

# plotHistCharge(ITTotalChargeSLC,LCType="SLC")
# plotHistCharge(ITTotalChargeHLC,LCType="HLC")





# import tables
# import pandas



# CORSIKA_ID = "DAT059871"
# outputDir = "/home/enpaudel/icecube/triggerStudy/simFiles/"


# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetTrigCount.hdf5", help='Input files after counting triggers')
# args = parser.parse_args()



# f = tables.open_file("/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/DAT059871GenDetTrigCount.hdf5")
# print(f.root)
# print("trig count", f.root)
# # print("trig count", trigCount.root.ITSMTTriggered.cols.Event[:])
# # trigCount



def plotHistCharge2(charges):
	chargeName = [k for k,v in locals().iteritems() if v == charges][0]
	if "SLC" in chargeName:
		LCType="SLC"
	if "HLC" in chargeName:
		LCType="HLC"
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	chargeMean = np.mean(charges)
	chargeSD = np.std(charges)
	minLim = min(charges)
	maxLim = max(charges)
	n_bins = np.linspace(minLim,maxLim,20)
	# n_bins = np.linspace(0,90,9)
	ax.hist(charges,bins=n_bins,histtype="step",lw=2.5,edgecolor='#2ca02c',label="hist",alpha=1)
	ax.axvline(x=chargeMean, ymin=0, ymax=1,c='#2ca02c',ls="-.",label=r"$\bar{q}$")
	ax.axvspan(chargeMean - chargeSD,chargeMean + chargeSD,ymin=0, ymax=1,color='#2ca02c',ls="-",label=r"$\sigma$",alpha=0.1)
	# ax.hist(thetaList,bins=n_bins,histtype="step",lw=2.5,fc='orange',label="ratio: 12 - 75",alpha=0.75)
	ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
	ax.tick_params(axis='both',which='major', length=8,width=1.2)
	ax.tick_params(axis='both',which='minor', length=4)
	# ax.set_yscale('log')
	ax.grid(True)
	# ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
	# ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
	ax.set_xlabel(r"q [vem]", fontsize=28)
	ax.set_ylabel(r"$count$", fontsize=22)
	ax.set_ylim(1,None)
	ax.legend(fontsize=20)
	plt.savefig(plotFolder+"/histCharge"+str(LCType)+".pdf",transparent=False,bbox_inches='tight')
	plt.close()