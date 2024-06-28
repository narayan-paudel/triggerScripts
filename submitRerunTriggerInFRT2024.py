#!/usr/bin/env python3

import os
import glob
import subprocess

import json
import datetime as dt
from astropy.time import Time


from dateutil import parser

import numpy as np
import pandas as pd

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"

year = "2024"
submitFileName = ABS_PATH_HERE+"calcRate2024FRT.sub"

def makeSubFile(fileList):
  # print("writing submit file for corsika",corsikaFile)
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/rerunTriggerInFRT2024.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/FRT2024{}.log\n".format(year))
  submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/FRT2024{}.out\n".format(year))
  submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/FRT2024{}.err\n".format(year))
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
  submitFile.write("arguments = ")
  for ifile in fileList:
    submitFile.write(ifile + " ")
  submitFile.write("\n")
  submitFile.write("queue 1\n")
  submitFile.close()

def submitToCondorFile(fileList,run):
  makeSubFile(fileList)
  subprocess.call(["condor_submit calcRate2024FRT.sub -batch-name {0}".format(run)], shell=True)
  subprocess.call(["rm calcRate2024FRT.sub"], shell=True)




##############for all 2023 runs season#############################
###################################################################
# f = open('/home/enpaudel/icecube/triggerStudy/GRL2023/GRL2023.json')
# data = json.load(f)
# for i in range(len(data["runs"])):
#   if data["runs"][i]["good_it"]:
#     run = data["runs"][i]["run"]
#     t_start = pd.to_datetime(data["runs"][i]["good_tstart"], format='%Y-%m-%d %H:%M:%S.%f')
#     t_stop = pd.to_datetime(data["runs"][i]["good_tstop"], format='%Y-%m-%d %H:%M:%S.%f')
#     # datetime.strptime(datetime_str, '%m/%d/%y %H:%M:%S')
#     # print(data["runs"][i]["good_tstart"])
#     # t_start = Time(data["runs"][i]["good_tstart"])
#     # t_stop = Time(data["runs"][i]["good_tstop"])
#     # # t_start = dataclasses.I3Time(data["runs"][i]["good_tstart"])
#     # # t_stop = dataclasses.I3Time(data["runs"][i]["good_tstop"])
#     # t_start = parser.parse(data["runs"][i]["good_tstart"])
#     # t_stop = parser.parse(data["runs"][i]["good_tstop"])
#     # t_start = dt.datetime.strptime(data["runs"][i]["good_tstart"],'%Y-%m-%d %H:%M:%S.%f')
#     # t_stop = dt.datetime.strptime(data["runs"][i]["good_tstop"],'%Y-%m-%d %H:%M:%S.%f')
#     # print("info",run, t_stop.to_value('mjd', 'long'), t_start.to_value('mjd', 'long'))
#     # print("info",run, t_stop.to_value('mjd', 'long')-t_start.to_value('mjd', 'long'))
#     # print("info",run, t_stop.sec, t_start.sec)
#     print("info",run, (t_stop-t_start).total_seconds())

###################################################################
######################for rusn in GRL2023##########################
GRL = "/data/exp/IceCube/2023/filtered/level2/"
GRL2023 = GRL + "IC86_2023_GoodRunInfo.txt"
GRL2023_df = pd.read_csv(GRL2023,sep="\\s+",header=0,escapechar='#')
GRL_IT_List = np.asarray(GRL2023_df.iloc[:,0][GRL2023_df.iloc[:,2]==1])
liveTime = np.asarray(GRL2023_df.iloc[:,3][GRL2023_df.iloc[:,2]==1])
runList = [int(irun) for irun in GRL_IT_List]
for irun in runList:
  if irun < 138808:
    fileList = sorted(glob.glob("/data/exp/IceCube/2023/filtered/PFFilt/1[12]??/*{}*".format(irun)))
  elif irun >= 138808:
    fileList = sorted(glob.glob("/data/exp/IceCube/2024/filtered/PFFilt/0[1234]??/*{}*".format(irun)))
  submitToCondorFile(fileList,irun)


# print(fileList)




# print(runList)
# data2024 = "/data/exp/IceCube/2024/filtered/PFFilt/"







############################rest###################################
###################################################################

# runs = np.linspace(138858,138882,25)
# # runs = [138867]
# runs = [int(irun) for irun in runs]
# data2024 = "/data/exp/IceCube/2024/filtered/PFFilt/"
# for irun in runs:
#   fileList = sorted(glob.glob(data2024+"01??/*{}*".format(irun)))
#   submitToCondorFile(fileList,irun)


#for Decemeber set of runs
# runs_highRate = np.linspace(138753,138774,22)
# runs_highRate = [int(irun) for irun in runs_highRate]
# outlierEvents = [138764,138768]
# for iout in outlierEvents:
#   runs_highRate.remove(iout)
# print(runs_highRate)
# # runs_highRate = [138765] #6GB 
# # runs_highRate = [138771] #9GB
# data2023 = "/data/exp/IceCube/2023/filtered/PFFilt/"
# for irun in runs_highRate:
#   fileList = sorted(glob.glob(data2023+"12??/*{}*".format(irun)))
#   submitToCondorFile(fileList,irun)

# #for 2023 January runs ( January 13 and January 21)
# runs = np.linspace(137540,137570,31)
# runs = [int(irun) for irun in runs]
# outlierEvents = [137562,137555,137556,137557]
# for iout in outlierEvents:
#   runs.remove(iout)
# print(runs)
# data2023 = "/data/exp/IceCube/2023/filtered/PFFilt/"
# for irun in runs:
#   fileList = sorted(glob.glob(data2023+"01??/*{}*".format(irun)))
#   submitToCondorFile(fileList,irun)



# #for 2022 January runs ( January 13 and January 21)
# runs = np.linspace(136163,136190,28)
# runs = [int(irun) for irun in runs]
# outlierEvents = [136167]
# for iout in outlierEvents:
#   runs.remove(iout)
# print(runs)
# data2022 = "/data/exp/IceCube/2022/filtered/PFFilt/"
# for irun in runs:
#   fileList = sorted(glob.glob(data2022+"01??/*{}*".format(irun)))
#   submitToCondorFile(fileList,irun)


#for 2021 January runs ( January 13 and January 21)
# runs = np.linspace(134888,134915,28)
# runs = [int(irun) for irun in runs]
# outlierEvents = [134892]
# for iout in outlierEvents:
#   runs.remove(iout)
# print(runs)
# runs = [134893] #9 GB
# data2021 = "/data/exp/IceCube/2021/filtered/PFFilt/"
# for irun in runs:
#   fileList = sorted(glob.glob(data2021+"01??/*{}*".format(irun)))
#   submitToCondorFile(fileList,irun)


