#!/usr/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
basePathProton = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/"
basePathIron = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10889/"
basePathHelium = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/11663/"
basePathOxygen = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/12605/"
print("abs path",ABS_PATH_HERE)
#no of files 10410 = 59997, else 60000
energyList = ['5.0', '5.1', '5.2', '5.3', '5.4', '5.5', '5.6', '5.7', '5.8', '5.9', '6.0', '6.1', '6.2', '6.3', '6.4', '6.5', '6.6', '6.7', '6.8', '6.9', '7.0', '7.1', '7.2', '7.3', '7.4', '7.5', '7.6', '7.7', '7.8', '7.9']
# energyList = ['7.3']

simFiles = "/home/enpaudel/icecube/triggerStudy/simFiles/"
#need corresponding change of dataSetUnique folder in runDetectorL2Cluster, addUniqueRunID.py and DetectorSim.py
# dataSetUnique = "dataSetUnique/"
# dataSetUnique = "dataSetUniqueFRT/"
dataSetUnique = "dataSetUniqueWFRT/"

print("energy list",energyList)
submitFileName = ABS_PATH_HERE+"tempSubmitL2.sub"

def makeSubFile(corsikaFile,primary):
	# print("writing submit file for corsika",corsikaFile)
	corsikaID = str(corsikaFile).split("/")[-1]
	energyID = float(str(corsikaFile).split("/")[-2])
	# print("energyID",energyID)
	corsikaID = corsikaID.split(".")[0]
	submitFile = open(submitFileName,"w")
	submitFile.write("########################################\n")
	submitFile.write("## submit description file\n")
	submitFile.write("########################################\n\n")
	submitFile.write("Universe   = vanilla\n")
	submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/runDetectorL2Cluster.sh\n")
	submitFile.write("Log        = /scratch/enpaudel/log/trigStudyL2{0}{1}.log\n".format(primary,corsikaID))
	submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/trigStudyL2{0}{1}.out\n".format(primary,corsikaID))
	submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/trigStudyL2{0}{1}.err\n".format(primary,corsikaID))
	submitFile.write("request_cpus = 1\n")
	submitFile.write("request_memory = 3GB\n")
	submitFile.write("request_disk = 1GB\n")
	submitFile.write("#request_gpus = 1\n")
	submitFile.write("#should_transfer_files   = IF_NEEDED\n")
	submitFile.write("#when_to_transfer_output = ON_EXIT\n")
	submitFile.write("#notification = Complete\n")
	submitFile.write("#notify_user = <email-address>\n")
	submitFile.write("#priority = <integer>\n")
	submitFile.write("##long job\n")
	priority=100000
	# submitFile.write("+AccountingGroup=\"2_week.$ENV(USER)\" \n\n")
	# if energyID > 7.7:
	# 	# priority=0
	# 	submitFile.write("+AccountingGroup=\"2_week.$ENV(USER)\" \n\n")		
	# 	# submitFile.write('+AccountingGroup=\"long.$ENV(USER)\" \n')
	# 	# submitFile.write('+AccountingGroup = "long.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
	# 	# submitFile.write('+AccountingGroup = "1_week.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
	# 	# submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
	# elif energyID > 7.2:
	# 	submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
	# elif energyID > 7.0:
	# 	submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
	submitFile.write("priority = {}\n".format(priority))
	submitFile.write("#set arguments to executable\n")
	submitFile.write("arguments = {0}\n\n".format(corsikaFile))
	submitFile.write("queue 1\n")
	submitFile.close()

def getCORSIKALists(basePath,energyDir):
	corsikaList=sorted(glob.glob(basePath+str(energyDir)+"/*.bz2"))
	return corsikaList

def getCorsikaFiles(basePath,energyList):
	# energyList=sorted([f.name for f in os.scandir(basePath) if f.is_dir()])
	print("energyDir",energyList)
	corsikaFiles=[]
	for ienergy in energyList:
		print("ienergy",ienergy)
		corsikaList=getCORSIKALists(basePath,ienergy)
		# print("corsikaList",corsikaList)
		# corsikaFiles.append(corsikaList[:11])
		corsikaFiles += corsikaList[0:200]
		# corsikaFiles += corsikaList
	# print("corsika files",corsikaFiles)
	return corsikaFiles

def submitToCondorFile(corsikaFile,primary):
	makeSubFile(corsikaFile,primary)
	corsikaID = str(corsikaFile).split("/")[-1]
	corsikaID = str(corsikaID).split(".")[0]
	simulatedGenFile = simFiles+"dataSetGen/"+primary+corsikaID+"Gen.i3.bz2"
	simulatedUniqueFile = simFiles+dataSetUnique+primary+corsikaID+"GenDetFiltProcUnique.i3.gz"
	# print("check files",simulatedGenFile,simulatedUniqueFile,os.path.exists(simulatedGenFile), not simulatedUniqueFile,os.path.exists(simulatedGenFile) and not simulatedUniqueFile)
	if os.path.exists(simulatedGenFile) and not os.path.exists(simulatedUniqueFile):
		corsikaID = int(''.join(i for i in corsikaID if i.isdigit()))
		batName = "{0}--{1}".format(primary,corsikaID)
		subprocess.call(["condor_submit tempSubmitL2.sub -batch-name {0}--{1}".format(primary,corsikaID)], shell=True)
		print("batch name",batName)
		print('primary',primary,corsikaID)
	subprocess.call(["rm tempSubmitL2.sub"], shell=True)
def submitToCondor(corsikaFiles,primary):	
	for ifiles in corsikaFiles:
		submitToCondorFile(ifiles,primary)

# # print("submitting proton files: ", basePathProton)
# print("submitting Iron files from: ", basePathIron)
# ironCorsikaFiles = getCorsikaFiles(basePathIron,energyList)
# submitToCondor(ironCorsikaFiles,"Fe")

# protonCorsikaFiles = getCorsikaFiles(basePathProton,energyList)
# submitToCondor(protonCorsikaFiles,"p")

# heliumCorsikaFiles = getCorsikaFiles(basePathHelium,energyList)
# submitToCondor(heliumCorsikaFiles,"He")

oxygenCorsikaFiles = getCorsikaFiles(basePathOxygen,energyList)
submitToCondor(oxygenCorsikaFiles,"O")
