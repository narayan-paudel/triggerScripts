#!/usr/bin/env python3

import os
import glob
import subprocess

import numpy as np
from inclinedTriggerTools import nCorFiles

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"

# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/"
basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/"
# basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTanks/"

# hdf5NullListP = sorted(glob.glob(basePath+"p*Clean*.hdf5"))
# hdf5NullListHe = sorted(glob.glob(basePath+"He*Clean*.hdf5"))
# hdf5NullListO = sorted(glob.glob(basePath+"O*Clean*.hdf5"))
# hdf5NullListFe = sorted(glob.glob(basePath+"Fe*Clean*.hdf5"))
# hdf5NullList = np.concatenate((hdf5NullListP,hdf5NullListHe,hdf5NullListO,hdf5NullListFe))

# print("file len H He O Fe",nCorFiles(hdf5NullListP),nCorFiles(hdf5NullListHe),nCorFiles(hdf5NullListO),nCorFiles(hdf5NullListFe))

hdf5List = sorted(glob.glob(basePath+"*Clean*.hdf5"))
print(hdf5List[:1])

submitFileName = ABS_PATH_HERE+"tempSubmitPickle.sub"

def makeSubFile(ihdf):
  # print("writing submit file for corsika",corsikaFile)
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/pickleEvtList.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/pickle.log\n")
  submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/pickle.out\n")
  submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/pickle.err\n")
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
  #   # priority=0
  #   submitFile.write("+AccountingGroup=\"2_week.$ENV(USER)\" \n\n")   
  #   # submitFile.write('+AccountingGroup=\"long.$ENV(USER)\" \n')
  #   # submitFile.write('+AccountingGroup = "long.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
  #   # submitFile.write('+AccountingGroup = "1_week.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
  #   # submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  # elif energyID > 7.2:
  #   submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  # elif energyID > 7.0:
  #   submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  submitFile.write("priority = {}\n".format(priority))
  submitFile.write("#set arguments to executable\n")
  submitFile.write("arguments = {0}\n\n".format(ihdf))
  submitFile.write("queue 1\n")
  submitFile.close()

def submitToCondor(ihdf):
  makeSubFile(ihdf)
  subprocess.call(["condor_submit tempSubmitPickle.sub"], shell=True)
  subprocess.call(["rm tempSubmitPickle.sub"], shell=True)


for ihdf in hdf5List:
  print(ihdf)
  submitToCondor(ihdf)



