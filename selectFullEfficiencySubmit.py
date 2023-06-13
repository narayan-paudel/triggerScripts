#!/bin/env python3

import os
import glob
import subprocess

import glob

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"

i3Files = sorted(glob.glob("/data/sim/IceCubeUpgrade/CosmicRay/Narayan/dataSetCleanTanks/*.i3.gz"))

print("input files",i3Files[:2])

submitFileName = ABS_PATH_HERE+"tempSubmitFullEff.sub"

def makeSubFile(i3File):
  # print("writing submit file for corsika",corsikaFile)
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/selectFullEfficiency.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/fullEff.log\n")
  submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/fullEff.out\n")
  submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/fullEff.err\n")
  submitFile.write("request_cpus = 1\n")
  submitFile.write("request_memory = 3GB\n")
  submitFile.write("request_disk = 1GB\n")
  submitFile.write("#request_gpus = 1\n")
  submitFile.write("#should_transfer_files   = IF_NEEDED\n")
  submitFile.write("#when_to_transfer_output = ON_EXIT\n")
  submitFile.write("#notification = Complete\n")
  submitFile.write("#notify_user = <email-address>\n")
  submitFile.write("#priority = <integer>\n")
  priority = 1000000000
  submitFile.write("##long job\n")
  # submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  submitFile.write("priority = {}\n".format(priority))
  submitFile.write("#set arguments to executable\n")
  submitFile.write("arguments = {0}\n\n".format(i3File))
  submitFile.write("queue 1\n")
  submitFile.close()

def submitToCondorFile(i3File):
  makeSubFile(i3File)
  subprocess.call(["condor_submit tempSubmitFullEff.sub"], shell=True)
  subprocess.call(["rm tempSubmitFullEff.sub"], shell=True)

def submitToCondor(fileList): 
  for ifiles in fileList:
    submitToCondorFile(ifiles)
submitToCondor(i3Files)
# submitToCondor(inFileSMT)
