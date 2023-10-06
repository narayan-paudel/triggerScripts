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

############################ if need to check GRL#
GRL = "/data/ana/CosmicRay/IceTop_GRL/"
GRL2020 = GRL + "/IC86_2020_GoodRunInfo_4IceTop_level2.txt"
GRL2020_df = pd.read_csv(GRL2020,sep="\t", header=0,escapechar='#')
GRL_IT_List = np.asarray(GRL2020_df["OutDir"][GRL2020_df["Good_it"]==1])
GRL_RunNum = np.asarray(GRL2020_df["RunNum"][GRL2020_df["Good_it"]==1])
GRL_LiveTime = np.asarray(GRL2020_df["LiveTime"][GRL2020_df["Good_it"]==1])
# print("GRL_RunNum",GRL_RunNum)

Level3_2021_path = "/data/ana/CosmicRay/IceTop_level3/exp/IC86.2020/2021/"
RunListFiles = sorted(glob.glob(Level3_2021_path+"03??/*/*.i3.*"))
RunList = sorted(glob.glob(Level3_2021_path+"03??/*"))
RunList = [int(irun.split("/")[-1].replace("Run00","")) for irun in RunList]
RunList = [irun for irun in RunList if irun in GRL_RunNum]
LiveTimeList = [itime for irun,itime in zip(GRL_RunNum,GRL_LiveTime) if irun in RunList]
totalTime = sum(LiveTimeList)

# print("runlist",len(RunList),RunList)
# print("timelist",len(LiveTimeList),sum(LiveTimeList)/3600.0,LiveTimeList)
###############################################
#for simulation
Level3_2021_sim_path = "/data/sim/IceTop/2023/generated/untriggered/level3"
SimListFiles = sorted(glob.glob(Level3_2021_sim_path+"/*.i3.gz"))
print(SimListFiles)
hdf5FileSimu = "/home/enpaudel/icecube/triggerStudy/plots/hdf5/verificationh5Simu.h"
hdf5FileData = "/home/enpaudel/icecube/triggerStudy/plots/hdf5/verificationh5data.h"


###############################################

RunListFiles = sorted(glob.glob(Level3_2021_path+"03??/*/*.i3.*"))
liveTime = glob.glob(Level3_2021_path+"03??/*/*.pickle")[0]
# with (open(liveTime, "rb")) as timefile:
#   fobj = pickle.load(timefile)
# print("fobj",fobj)

# RunList = sorted(glob.glob(Level3_2021_path+"03??/*"))
# print("runlist",len(RunList),RunList)








submitFileName = ABS_PATH_HERE+"tempSubmitL3Verify.sub"

def makeSubFile(fileList,isSim):
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/L3VerificationPlots.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/L3Verify$(Process).log\n")
  submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/L3Verify$(Process).out\n")
  submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/L3Verify$(Process).err\n")
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
  # submitFile.write("-i ")
  for ifile in fileList:
    submitFile.write(ifile + " ")
  if isSim:
    submitFile.write("--isSim")
  submitFile.write("\n")
  submitFile.write("queue 1\n")
  submitFile.close()



def submitToCondorFile(fileList,isSim):
  makeSubFile(fileList,isSim)
  subprocess.call(["condor_submit tempSubmitL3Verify.sub"], shell=True)
  subprocess.call(["rm tempSubmitL3Verify.sub"], shell=True)

def submitToCondor(fileList,chunk,isSim):
  fileChunks = [fileList[i:i + chunk] for i in range(0, len(fileList), chunk)]
  for ifileList in fileChunks:
    print("submitting to condor",ifileList[0])
    submitToCondorFile(ifileList,isSim)

submitToCondor(SimListFiles,5,isSim=True) #172 total
submitToCondor(RunListFiles,10,isSim=False) #2250 total
