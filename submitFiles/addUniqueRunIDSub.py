#!/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
print("abs path",ABS_PATH_HERE)
# inputPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSet/"
inputPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSet/"
inputList = sorted(glob.glob(inputPath+"*DAT*GenDetFiltProc.i3.bz2"))

submitFileName = ABS_PATH_HERE+"tempSubmitUnique.sub"


def makeSubFile(fileList):
	submitFile = open(submitFileName,"w")
	submitFile.write("########################################\n")
	submitFile.write("## submit description file\n")
	submitFile.write("########################################\n\n")
	submitFile.write("Universe   = vanilla\n")
	submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/addUniqueRunID.sh\n")
	submitFile.write("Log        = /scratch/enpaudel/log/Uniq$(Process).log\n")
	submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/Uniq$(Process).out\n")
	submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/Uniq$(Process).err\n")
	submitFile.write("request_cpus = 1\n")
	submitFile.write("request_memory = 2GB\n")
	submitFile.write("request_disk = 1GB\n")
	submitFile.write("#request_gpus = 1\n")
	submitFile.write("#should_transfer_files   = IF_NEEDED\n")
	submitFile.write("#when_to_transfer_output = ON_EXIT\n")
	submitFile.write("#notification = Complete\n")
	submitFile.write("#notify_user = <email-address>\n")
	submitFile.write("#priority = <integer>\n")
	submitFile.write("##long job\n")
	priority=1100
	# submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
	submitFile.write("priority = {}\n".format(priority))
	submitFile.write("#set arguments to executable\n")
	submitFile.write("arguments = ")
	for ifile in fileList:
		submitFile.write(ifile + " ")
	submitFile.write("\n")
	submitFile.write("queue 1\n")
	submitFile.close()

def submitToCondorFile(fileList):
	makeSubFile(fileList)
	subprocess.call(["condor_submit tempSubmitUnique.sub"], shell=True)
	# subprocess.call(["rm tempSubmit.sub"], shell=True)

def submitToCondor(fileList,chunk):
	fileChunks = [fileList[i:i + chunk] for i in range(0, len(fileList), chunk)]
	for ifileList in fileChunks:
		submitToCondorFile(ifileList)

submitToCondor(inputList,100)
