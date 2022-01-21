#!/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE=str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE="./"
ABS_PATH_HERE+="/"
basePath = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/"
print("abs path",ABS_PATH_HERE)
print("basePath",basePath)
energyList=sorted([f.name for f in os.scandir(basePath) if f.is_dir()])

def getCORSIKALists(energyDir):
	corsikaList=sorted(glob.glob(basePath+str(energyDir)+"/*.bz2"))
	return corsikaList

print("energyDir",energyList)
corsikaFiles=[]
for ienergy in energyList:
	corsikaList=getCORSIKALists(ienergy)
	corsikaFiles.append(corsikaList[0])
print("corsika files",corsikaFiles)

def submitToCondor(corsikaFile):
	subprocess.call(["sed 's|$(Item)|{}|' submitRunDetector.submit > tempSub.sub".format(corsikaFile)],shell=True)
	# subprocess.call(["condor_submit tempSub.sub -batch-name {0}".format(str(corsikaFile)).split("/")[-1]], shell=True)
	subprocess.call(["condor_submit tempSub.sub -batch-name {0}".format(str(corsikaFile).split("/")[-1])], shell=True)
	subprocess.call(["rm tempSub.sub"], shell=True)

for ifiles in corsikaFiles:
	submitToCondor(ifiles)