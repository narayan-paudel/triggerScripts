#!/bin/env python

import os
import glob
import subprocess

import pandas as pd
import numpy as np

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
print("abs path",ABS_PATH_HERE)

# inputFiles = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/"
inputFiles = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/"
pFiles = sorted(glob.glob(inputFiles+"pDAT*GenDetFiltProcUniqueCleanVEMEvts.i3.gz"))
HeFiles = sorted(glob.glob(inputFiles+"HeDAT*GenDetFiltProcUniqueCleanVEMEvts.i3.gz"))
OFiles = sorted(glob.glob(inputFiles+"ODAT*GenDetFiltProcUniqueCleanVEMEvts.i3.gz"))
FeFiles = sorted(glob.glob(inputFiles+"FeDAT*GenDetFiltProcUniqueCleanVEMEvts.i3.gz"))
submitFileName = ABS_PATH_HERE+"tempSubmitReco125.sub"

def makeSubFile(fileList):
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/reconstructIceTopS125.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/S125$(Process).log\n")
  submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/S125$(Process).out\n")
  submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/S125$(Process).err\n")
  submitFile.write("request_cpus = 1\n")
  submitFile.write("request_memory = 2GB\n")
  submitFile.write("request_disk = 2GB\n")
  submitFile.write("#request_gpus = 1\n")
  submitFile.write("#should_transfer_files   = IF_NEEDED\n")
  submitFile.write("#when_to_transfer_output = ON_EXIT\n")
  submitFile.write("#notification = Complete\n")
  submitFile.write("#notify_user = <email-address>\n")
  submitFile.write("#priority = <integer>\n")
  submitFile.write("##long job\n")
  priority=11000000000
  # submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  submitFile.write("priority = {}\n".format(priority))
  # submitFile.write("#set arguments to executable\n")
  submitFile.write("arguments = ")
  for ifile in fileList:
    submitFile.write(ifile + " ")
  submitFile.write("\n")
  submitFile.write("queue 1\n")
  submitFile.close()


def submitToCondorFile(fileList):
  makeSubFile(fileList)
  subprocess.call(["condor_submit tempSubmitReco125.sub"], shell=True)
  # subprocess.call(["rm tempSubmitS125.sub"], shell=True)

submitToCondorFile(pFiles)
submitToCondorFile(HeFiles)
submitToCondorFile(OFiles)
submitToCondorFile(FeFiles)