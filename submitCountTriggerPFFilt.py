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
# print("ITList",GRL_IT_List)
liveTime = np.asarray(GRL2023_df.iloc[:,3][GRL2023_df.iloc[:,2]==1])
# print(GRL2023_df.head())
runList = [int(irun) for irun in GRL_IT_List]
runList = [138615,138616,138617,138618,138619,138620,138621,138622,
138623,138624, 138626,138627,138628,138630,138631,138632,138633,
138634,138635,138636, 138637,138638,138647,138648,138649,138650,
138651,138652,138653,138654, 138655,138656,138657,138658,138660,
138661,138662,138663,138664,138665, 138666,138667,138668,138669,
138670,138671,138674,138675,138677,138678, 138680,138681,138682,
138683,138684,138685,138686,138687,138688,138689, 138690,138691,
138692,138693,138694,138695,138696,138697,138698,138699, 138700,
138701,138702,138703,138704,138706,138707,138708,138709,138710,
 138711,138729,138734,138735,138739,138740,138741,138742,138743,
 138751, 138752,138753,138754,138755,138756,138757,138758,138759,
 138760,138761, 138762,138763,138764,138765,138766,138767,138768,
 138769,138770,138771, 138772,138773,138774,138782,138783,138784,
 138786,138787,138788,138793, 138794,138801,138802,138803,138804,
 138805,138806,138807,138808,138809, 138810,138811,138812,138813,
 138814,138815,138816,138817,138818,138819, 138820,138821,138822,
 138823,138824,138826,138827,138828,138829,138830, 138831,138832,
 138833,138834,138835,138836,138837,138847,138848,138849, 138850,
 138851,138852,138853,138854,138855,138856,138857,138858,138859,
  138860,138861,138862,138863,138864,138865,138866,138867,138868,138869,
   138870,138871,138872,138873,138874,138875,138876,138877,138878,138879, 138880,138881,138882]
runList = [138615,138883]
runList = [139196,139197,139198] #for forbush decrease
chunk = 1
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