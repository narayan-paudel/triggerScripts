#!/bin/env python3

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
# energyList = ['5.0']#, '5.1', '5.2', '5.3', '5.4', '5.5', '5.6', '5.7', '5.8', '5.9', '6.0', '6.1', '6.2', '6.3', '6.4', '6.5', '6.6', '6.7', '6.8', '6.9', '7.0', '7.1', '7.2', '7.3', '7.4', '7.5', '7.6', '7.7', '7.8', '7.9']
energyList = ['5.0', '5.1', '5.2', '5.3', '5.4', '5.5', '5.6', '5.7', '5.8', '5.9', '6.0', '6.1', '6.2', '6.3', '6.4', '6.5', '6.6', '6.7', '6.8', '6.9', '7.0', '7.1', '7.2', '7.3', '7.4', '7.5', '7.6', '7.7', '7.8', '7.9']
# energyList = ['7.0', '7.1', '7.2', '7.3', '7.4', '7.5', '7.6', '7.7', '7.8', '7.9']

simFiles = "/home/enpaudel/icecube/triggerStudy/simFiles/"

print("energy list",energyList)
submitFileName = ABS_PATH_HERE+"tempSubmit.sub"


def makeSubFile(corsikaFile):
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
	submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/runDetectorCluster.sh\n")
	submitFile.write("Log        = /scratch/enpaudel/log/trigStudy{0}.log\n".format(corsikaID))
	submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/trigStudy{0}.out\n".format(corsikaID))
	submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/trigStudy{0}.err\n".format(corsikaID))
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
	# priority=4
	if energyID > 7.7:
		# priority=0
		submitFile.write("+AccountingGroup=\"2_week.$ENV(USER)\" \n\n")		
		# submitFile.write('+AccountingGroup=\"long.$ENV(USER)\" \n')
		# submitFile.write('+AccountingGroup = "long.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
		# submitFile.write('+AccountingGroup = "1_week.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
		# submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
	elif energyID > 7.2:
		submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
	# submitFile.write("priority = {}\n".format(priority))
	submitFile.write("#set arguments to executable\n")
	submitFile.write("arguments = {0}\n\n".format(corsikaFile))
	submitFile.write("queue 1\n")
	submitFile.close()

def runningJobs(filePath):
	with open(filePath,"r") as f:
		next(f)
		next(f)
		lines = f.readlines()
	lines = [iline.split(" ")[1] for iline in lines if iline.split(" ")[0] == "enpaudel" ]
	print("lines",lines)
	# lines = [iline.split("/")[-1] for iline in lines ]
	# lines = [iline.split(".")[0] for iline in lines ]
	# regex = re.compile(r"\d+")
	# lines = [int(x) for iline in lines for x in regex.findall(iline)]
	# print("lines",lines)
	return lines
runJobs = runningJobs("../runningJobs.txt")


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
		corsikaFiles += corsikaList[90:100]
		# corsikaFiles += corsikaList
	# print("corsika files",corsikaFiles)
	return corsikaFiles

def submitToCondorFile(corsikaFile,primary):
	makeSubFile(corsikaFile)
	corsikaID = str(corsikaFile).split("/")[-1]
	corsikaID = str(corsikaID).split(".")[0]
	simulatedFile = simFiles+"dataSetUnique/"+primary+corsikaID+"GenDetFiltProcUnique.i3.gz"
	# simulatedFile = simFiles+"dataSetUnique/"+primary+corsikaID+"GenDetFiltProcUnique.i3.gz"
	# print("simulated file",simulatedFile)
	if not os.path.exists(simulatedFile):
		corsikaID = int(''.join(i for i in corsikaID if i.isdigit()))
		batName = "{0}-{1}".format(primary,corsikaID)
		if batName not in runJobs:
			subprocess.call(["condor_submit tempSubmit.sub -batch-name {0}-{1}".format(primary,corsikaID)], shell=True)
			print('primary',primary,corsikaID)
	# subprocess.call(["rm tempSubmit.sub"], shell=True)
def submitToCondor(corsikaFiles,primary):	
	for ifiles in corsikaFiles:
		submitToCondorFile(ifiles,primary)

# print("submitting proton files: ", basePathProton)
protonCorsikaFiles = getCorsikaFiles(basePathProton,energyList)
submitToCondor(protonCorsikaFiles,"p")

heliumCorsikaFiles = getCorsikaFiles(basePathHelium,energyList)
submitToCondor(heliumCorsikaFiles,"He")

oxygenCorsikaFiles = getCorsikaFiles(basePathOxygen,energyList)
submitToCondor(oxygenCorsikaFiles,"O")

# print("submitting Iron files from: ", basePathIron)
ironCorsikaFiles = getCorsikaFiles(basePathIron,energyList)
submitToCondor(ironCorsikaFiles,"Fe")