#!/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
basePath = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/"
print("abs path",ABS_PATH_HERE)
print("basePath",basePath)
energyList=sorted([f.name for f in os.scandir(basePath) if f.is_dir()])
print("energy list",energyList)
# energyList=energyList[:1]
# energyList = ['5.0', '5.1', '5.2', '5.3', '5.4', '5.5', '5.6', '5.7', '5.8', '5.9', '6.0', '6.1', '6.2', '6.3', '6.4', '6.5', '6.6', '6.7', '6.8', '6.9', '7.0', '7.1', '7.2', '7.3', '7.4', '7.5', '7.6', '7.7', '7.8', '7.9']
# energyList = ['5.0']#, '5.1', '5.2', '5.3', '5.4', '5.5', '5.6', '5.7', '5.8', '5.9', '6.0', '6.1', '6.2', '6.3', '6.4', '6.5', '6.6', '6.7', '6.8', '6.9', '7.0']
energyList = ['7.0', '7.1', '7.2', '7.3', '7.4', '7.5', '7.6', '7.7', '7.8', '7.9']
# energyList = ['5.0', '5.1', '5.2', '5.3', '5.4', '5.5', '5.6', '5.7', '5.8', '5.9', '6.0', '6.1', '6.2', '6.3', '6.4', '6.5', '6.6', '6.7', '6.8', '6.9']

print("energy list",energyList)
submitFileName = ABS_PATH_HERE+"tempSubmit.sub"


def makeSubFile(corsikaFile):
	print("writing submit file for corsika",corsikaFile)
	corsikaID = str(corsikaFile).split("/")[-1]
	energyID = float(str(corsikaFile).split("/")[-2])
	print("energyID",energyID)
	corsikaID = corsikaID.split(".")[0]
	submitFile = open(submitFileName,"w")
	submitFile.write("########################################\n")
	submitFile.write("## submit description file\n")
	submitFile.write("########################################\n\n")
	submitFile.write("Universe   = vanilla\n")
	submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/scripts/runDetectorCluster.sh\n")
	submitFile.write("Log        = /scratch/enpaudel/log/trigStudy{0}.log\n".format(corsikaID))
	submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/trigStudy{0}.out\n".format(corsikaID))
	submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/trigStudy{0}.err\n".format(corsikaID))
	submitFile.write("request_cpus = 1\n")
	submitFile.write("request_memory = 4GB\n")
	submitFile.write("request_disk = 1GB\n")
	submitFile.write("#request_gpus = 1\n")
	submitFile.write("#should_transfer_files   = IF_NEEDED\n")
	submitFile.write("#when_to_transfer_output = ON_EXIT\n")
	submitFile.write("#notification = Complete\n")
	submitFile.write("#notify_user = <email-address>\n")
	submitFile.write("#priority = <integer>\n")
	submitFile.write("##long job\n")
	if energyID > 6.9:
		# submitFile.write('#+AccountingGroup = "long.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
		# submitFile.write('+AccountingGroup = "1_week.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
		submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
	submitFile.write("#set arguments to executable\n")
	submitFile.write("arguments = {0}\n\n".format(corsikaFile))
	submitFile.write("queue 1\n")
	submitFile.close()


def getCORSIKALists(energyDir):
	corsikaList=sorted(glob.glob(basePath+str(energyDir)+"/*.bz2"))
	return corsikaList

print("energyDir",energyList)
corsikaFiles=[]
for ienergy in energyList:
	print("ienergy",ienergy)
	corsikaList=getCORSIKALists(ienergy)
	# print("corsikaList",corsikaList)
	# corsikaFiles.append(corsikaList[:11])
	corsikaFiles += corsikaList[:100]
# print("corsika files",corsikaFiles)

def submitToCondor(corsikaFile):
	makeSubFile(corsikaFile)
	subprocess.call(["condor_submit tempSubmit.sub -batch-name {0}".format(str(corsikaFile).split("/")[-1])], shell=True)
	# subprocess.call(["rm tempSubmit.sub"], shell=True)

for ifiles in corsikaFiles:
	submitToCondor(ifiles)