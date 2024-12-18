#!/usr/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
print("abs path",ABS_PATH_HERE)

# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('inputPath', type=str, default="/home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique/", help='Source folder for input files.')
# parser.add_argument('outputPath', type=str, default="/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/", help='Destination folder for output files.')
# args = parser.parse_args()
############################################################################
# inputPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique/"
# inputPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique1_6/"
inputPath = "/data/user/enpaudel/triggerStudy/simFiles/dataSetClean1_6/"
# inputPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique.FRT/"
# inputPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetUniqueSeedSame/"
OinputList = sorted(glob.glob(inputPath+"ODAT*GenDetFiltProcUnique*.i3.*"))
pinputList = sorted(glob.glob(inputPath+"pDAT*GenDetFiltProcUnique*.i3.*"))
HeinputList = sorted(glob.glob(inputPath+"HeDAT*GenDetFiltProcUnique*.i3.*"))
FeinputList = sorted(glob.glob(inputPath+"FeDAT*GenDetFiltProcUnique*.i3.*"))
###########################################################################
# inputPath = "/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/"
# OinputList = sorted(glob.glob(inputPath+"12631/"+"Level3_IC86.2012_12631_*.i3.gz"))[:3000]
# pinputList = sorted(glob.glob(inputPath+"12360/"+"Level3_IC86.2012_12360_*.i3.gz"))[:3000]
# HeinputList = sorted(glob.glob(inputPath+"12630/"+"Level3_IC86.2012_12630_*.i3.gz"))[:3000]
# FeinputList = sorted(glob.glob(inputPath+"12362/"+"Level3_IC86.2012_12362_*.i3.gz"))[:3000]
#############################################################################
print("inputList",pinputList)



submitFileName = ABS_PATH_HERE+"splitAndFilt.sub"


def makeSubFile(fileList):
	submitFile = open(submitFileName,"w")
	submitFile.write("########################################\n")
	submitFile.write("## submit description file\n")
	submitFile.write("########################################\n\n")
	submitFile.write("Universe   = vanilla\n")
	submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/runSplitAndFilter7HG.sh\n")
	submitFile.write("Log        = /scratch/enpaudel/log/splitAndFilt$(Process).log\n")
	submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/splitAndFilt$(Process).out\n")
	submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/splitAndFilt$(Process).err\n")
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
	priority=9900000
	# submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
	submitFile.write("priority = {}\n".format(priority))
	submitFile.write("#set arguments to executable\n")
	submitFile.write("arguments = ")
	for ifile in fileList:
		submitFile.write(ifile + " ")
	submitFile.write("\n")
	submitFile.write("queue 1\n")
	submitFile.close()

def submitToCondorFile(fileList,primary):
	makeSubFile(fileList)
	corsikaID = str(fileList[0]).split("/")[-1]
	corsikaID = corsikaID.split(".")[0]
	corsikaID = int(''.join(i for i in corsikaID if i.isdigit()))
	# print("file list",*fileList[:2])
	print("corsika id",corsikaID,primary)
	subprocess.call(["condor_submit splitAndFilt.sub -batch-name {0}---{1}".format(primary,corsikaID)], shell=True)
	subprocess.call(["rm splitAndFilt.sub"], shell=True)

def submitToCondor(fileList,chunk,primary):
	fileChunks = [fileList[i:i + chunk] for i in range(0, len(fileList), chunk)]
	for ifileList in fileChunks:
		print("submitting to condor",ifileList[0])
		submitToCondorFile(ifileList,primary)

submitToCondor(pinputList,1,"p")
submitToCondor(HeinputList,1,"He")
submitToCondor(OinputList,1,"O")
submitToCondor(FeinputList,1,"Fe")

