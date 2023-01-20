#!/usr/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"

submitFileName = ABS_PATH_HERE+"tempSubmitZenithReco.sub"

def makeSubFile():
  # print("writing submit file for corsika",corsikaFile)
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/applyZenithRecoFilter.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/trigStudyZenReco.log\n")
  submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/trigStudyZenReco.out\n")
  submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/trigStudyZenReco.err\n")
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
  submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  submitFile.write("priority = {}\n".format(priority))
  submitFile.write("#set arguments to executable\n")
  # submitFile.write("arguments = {0}\n\n".format(corsikaFile))
  submitFile.write("queue 1\n")
  submitFile.close()



def submitToCondor():
  makeSubFile()
  subprocess.call(["condor_submit tempSubmitZenithReco.sub"], shell=True)
  subprocess.call(["rm tempSubmitZenithReco.sub"], shell=True)


submitToCondor()
