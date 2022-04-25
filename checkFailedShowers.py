import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
from matplotlib.ticker import AutoMinorLocator,MultipleLocator
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap

import numpy as np
import subprocess
from customColors import qualitative_colors

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 20})

idList = np.linspace(1,3571,120)
idList = [int(n) for n in idList]
print("idList",idList)

logPath_base = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.0/"

logPathList = [logPath_base+"DAT{0:06d}.log".format(ID) for ID in idList]
print(logPathList[100])
plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"


def logZenith(logPath):
	''' extracts theta value from the given path'''
	theta = subprocess.check_output(['awk', '/PRIMARY ANGLES ARE:/ {print $6}', logPath])
	return np.rad2deg(float(theta.strip()))

def logPhi(logPath):
	''' extracts theta value from the given path'''
	phi = subprocess.check_output(['awk', '/PRIMARY ANGLES ARE:/ {print $10}', logPath])
	return np.rad2deg(float(phi.strip()))

def logCore(logPath):
	''' extracts theta value from the given path'''
	core = subprocess.check_output(['awk', '/PRIMARY ANGLES ARE:/ {print $10}', logPath])
	return np.rad2deg(float(phi.strip()))



thetaList = [logZenith(ipath) for ipath in logPathList]
phiList = [logPhi(ipath) for ipath in logPathList]

def plotThetaPhi(idList,thetaList,phiList):
	fig = plt.figure(figsize=(8,5))
	gs = gridspec.GridSpec(nrows=1,ncols=1)
	ax = fig.add_subplot(gs[0])
	ax.plot(idList,thetaList,"o-",lw=2.5,c=qualitative_colors(3)[1],label="theta",alpha=1)
	ax.plot(idList,phiList,"o-",lw=2.5,c=qualitative_colors(3)[2],label="phi",alpha=1)
	ax.axvspan(xmin=1171,xmax=1501,ymin=0,ymax=1,color="orange")
	# ax.plot(angleBins,totalEvts_list,"o-",c=next(colorsIter),label="total Evts",alpha=1)
	ax.set_xlabel(r"id", fontsize=20)
	ax.set_ylabel(r"angle [$^{\circ}$]", fontsize=20)
	ax.tick_params(axis='both',which='both',direction='in', labelsize=20)
	ax.tick_params(which='both', width=1.5)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	# ax.set_yscale('log')
	ax.grid(True,alpha=0.5)
	# ax.set_ylim(0,None)
	ax.legend(fontsize=14)
	# ax.legend(fontsize=14,ncol=2)
	# ax.yaxis.set_minor_locator(MultipleLocator(100))
	# ax.xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.savefig(plotFolder+"thetaPhi.pdf",transparent=False,bbox_inches='tight')
	plt.close()

plotThetaPhi(idList,thetaList,phiList)

