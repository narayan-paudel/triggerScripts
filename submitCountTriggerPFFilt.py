#!/usr/bin/env python3

import os
import glob
import subprocess

from numpy import linspace

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"

submitFileName = ABS_PATH_HERE+"countTrigger.sub"

def makeSubFile(GCD,fileList):
  # print("writing submit file for corsika",corsikaFile)
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/runCountTriggerPFFilt.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/countTrigPFFilt.log\n")
  submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/countTrigPFFilt.out\n")
  submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/countTrigPFFilt.err\n")
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
  priority=100000
  # submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  submitFile.write("priority = {}\n".format(priority))
  submitFile.write("#set arguments to executable\n")
  # submitFile.write("arguments = {0} {1}\n\n".format(GCD,subrun))
  submitFile.write("arguments = {} ".format(GCD))
  for ifile in fileList:
          submitFile.write(ifile + " ")
  submitFile.write("\n")
  submitFile.write("queue 1\n")
  submitFile.close()



def submitToCondor(GCD,inputList,run):
  print(inputList)
  makeSubFile(GCD,inputList)
  sub_run = inputList[0].split("_")[-1]
  sub_run = sub_run.split(".")[0]
  subprocess.call(["condor_submit countTrigger.sub -batch-name {0}-{1}".format(run,sub_run)], shell=True)
  subprocess.call(["rm countTrigger.sub"], shell=True)

pathTOPFFilt = "/data/exp/IceCube/2024/filtered/PFFilt/"
pathTOPFGCD = "/data/exp/IceCube/2024/internal-system/sps-gcd/"
pathTOGCD = "/data/sim/IceTop/2023/generated/untriggered/run2023/GCD/"

#####################use custom runs
# runList = linspace(138935,138952,18)
# runList = [int(irun) for irun in runList]
# chunk = 10
# for irun in runList:
#   # if not os.path.exists(GCD)
#   PFGCD = glob.glob(pathTOPFGCD+"*/PFGCD_Run00{}*".format(irun))[0]
#   subprocess.call(["tar -xvf {} -C  /data/sim/IceTop/2023/generated/untriggered/run2023/GCD".format(PFGCD)], shell=True)
#   GCD = glob.glob(pathTOGCD+"/PFGCD_Run00{}*.i3.gz".format(irun))[0]
#   inputList = sorted(glob.glob(pathTOPFFilt+"*/PFFilt_PhysicsFiltering_Run00{}*.tar.bz2".format(irun)))
#   fileChunks = [inputList[i:i + chunk] for i in range(0, len(inputList), chunk)]
#   for ifileChunks in fileChunks:
#     submitToCondor(GCD,ifileChunks,irun)
######################################using GRL##############
import pandas as pd
import numpy as np
GRL = "/data/exp/IceCube/2023/filtered/level2/"
GRL2023 = GRL + "IC86_2023_GoodRunInfo.txt"
GRL2023_df = pd.read_csv(GRL2023,sep="\\s+",header=0,escapechar='#')
# colNames = ["RunNum", "Good_i3", "Good_it", "LiveTime", "ActiveStrings", "ActiveDoms", "ActiveInIce", "OutDir                                                      ", "Comment(s)"]
# print("column names",GRL2023_df.columns)
# GRL2023_df[["RunNum", "Good_i3", "Good_it", "LiveTime", "ActiveStrings", "ActiveDoms", "ActiveInIce", "OutDir", "Comment(s)"]] = GRL2023_df["RunNum Good_i3 Good_it LiveTime ActiveStrings ActiveDoms ActiveInIce OutDir                                                       Comment(s) "].str.split(" ",expand=True)
GRL_IT_List = np.asarray(GRL2023_df.iloc[:,0][GRL2023_df.iloc[:,2]==1])
liveTime = np.asarray(GRL2023_df.iloc[:,3][GRL2023_df.iloc[:,2]==1])
# print(GRL2023_df.head())
runList = [int(irun) for irun in GRL_IT_List]
chunk = 100
for irun in runList:
  # if not os.path.exists(GCD)
  if irun >= 138808:
    pathTOPFFilt = "/data/exp/IceCube/2024/filtered/PFFilt/"
    pathTOPFGCD = "/data/exp/IceCube/2024/internal-system/sps-gcd/"
  elif irun < 138808:
    pathTOPFFilt = "/data/exp/IceCube/2023/filtered/PFFilt/"
    pathTOPFGCD = "/data/exp/IceCube/2023/internal-system/sps-gcd/"
  pathTOGCD = "/data/sim/IceTop/2023/generated/untriggered/run2023/GCD/"
  PFGCD = glob.glob(pathTOPFGCD+"*/PFGCD_Run00{}*".format(irun))[0]
  subprocess.call(["tar -xvf {} -C  /data/sim/IceTop/2023/generated/untriggered/run2023/GCD".format(PFGCD)], shell=True)
  GCD = glob.glob(pathTOGCD+"/PFGCD_Run00{}*.i3.gz".format(irun))[0]
  inputList = sorted(glob.glob(pathTOPFFilt+"*/PFFilt_PhysicsFiltering_Run00{}*.tar.bz2".format(irun)))
  fileChunks = [inputList[i:i + chunk] for i in range(0, len(inputList), chunk)]
  for ifileChunks in fileChunks:
    submitToCondor(GCD,ifileChunks,irun)