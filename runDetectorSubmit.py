#!/bin/env python3

import os
import glob
import subprocess

import numpy as np

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
basePathProton = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/"
basePathIron = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10889/"
basePathHelium = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/11663/"
basePathOxygen = "/data/sim/IceTop/2009/generated/CORSIKA-ice-top/12605/"
print("abs path",ABS_PATH_HERE)
#no of files 10410 = 59997, else 60000
energyList = ['5.0', '5.1', '5.2', '5.3', '5.4', '5.5', '5.6', '5.7', '5.8', '5.9', '6.0', '6.1', '6.2', '6.3', '6.4', '6.5', '6.6', '6.7', '6.8', '6.9', '7.0', '7.1', '7.2', '7.3', '7.4', '7.5', '7.6', '7.7', '7.8', '7.9']
# energyList = ['7.0', '7.1', '7.2', '7.3', '7.4', '7.5', '7.6', '7.7', '7.8', '7.9']

simFiles = "/home/enpaudel/icecube/triggerStudy/simFiles/"

print("energy list",energyList)
submitFileName = ABS_PATH_HERE+"tempSubmit.sub"


def makeSubFile(corsikaFile,primary):
  # print("writing submit file for corsika",corsikaFile)
  corsikaID = str(corsikaFile).split("/")[-1]
  energyID = float(str(corsikaFile).split("/")[-2])
  # print("energyID",energyID)
  corsikaID = corsikaID.split(".")[0]
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /home/enpaudel/icecube/triggerStudy/triggerScripts/runDetectorCluster.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/trigStudy{0}{1}.log\n".format(primary,corsikaID))
  submitFile.write("Output     = /data/user/enpaudel/triggerStudy/log/trigStudy{0}{1}.out\n".format(primary,corsikaID))
  submitFile.write("Error      = /data/user/enpaudel/triggerStudy/log/trigStudy{0}{1}.err\n".format(primary,corsikaID))
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
  # priority=100000
  # submitFile.write("+AccountingGroup=\"2_week.$ENV(USER)\" \n\n")
  # submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  if energyID > 7.7:
    # priority=0
    submitFile.write("+AccountingGroup=\"2_week.$ENV(USER)\" \n\n")   
    # submitFile.write('+AccountingGroup=\"long.$ENV(USER)\" \n')
    # submitFile.write('+AccountingGroup = "long.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
    # submitFile.write('+AccountingGroup = "1_week.$ENV(USER)" #other options 1_week, 2_week, instead of long\n')
    # submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  elif energyID > 7.5:
    submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  elif energyID > 7.0:
    submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  # submitFile.write("priority = {}\n".format(priority))
  submitFile.write("#set arguments to executable\n")
  submitFile.write("arguments = {0}\n\n".format(corsikaFile))
  submitFile.write("queue 1\n")
  submitFile.close()

def runningJobs(filePath):
  '''after running condor_q -run -batch > ../runningJobs.txt and condor_q -idle -batch > ../idleJobs.txt'''
  with open(filePath,"r") as f:
    next(f)
    next(f)
    lines = f.readlines()
  lines = [iline.split(" ")[1] for iline in lines if iline.split(" ")[0] == "enpaudel" ]
  # lines = [iline.split("/")[-1] for iline in lines ]
  # lines = [iline.split(".")[0] for iline in lines ]
  # regex = re.compile(r"\d+")
  # lines = [int(x) for iline in lines for x in regex.findall(iline)]
  # print("lines",lines)
  return lines
runJobs = runningJobs("../runningJobs.txt")

def idleJobs(filePath):
  with open(filePath,"r") as f:
    next(f)
    next(f)
    lines = f.readlines()
  lines = [iline.split(" ")[1].rstrip() for iline in lines if iline.split(" ")[0] == "enpaudel" ]
  # lines = [iline.split("/\")[-1] for iline in lines if iline.split(" ")[0] == "enpaudel" ]
  # lines = [iline.split("/")[-1] for iline in lines ]
  # lines = [iline.split(".")[0] for iline in lines ]
  # regex = re.compile(r"\d+")
  # lines = [int(x) for iline in lines for x in regex.findall(iline)]
  # print("lines",lines)
  return lines
idleLines = idleJobs("../idleJobs.txt")


queueJobs = runJobs+idleLines
queueJobs = []
# queueJobs = runJobs

jobs7 = []
jobs14 = ["Fe-867","Fe-897","Fe-807","Fe-778","Fe-808","He-777","He-897","He-778","Fe-1108",,"Fe-1138","He-808","He-838","O-57","Fe-87","He-87","O-117","Fe-237","O-237","Fe-267","Fe-297","Fe-627","Fe-687",
"Fe-747","Fe-628","Fe-658","Fe-1047","p-1047","p-1018","p-1050","He-928","O-1049","O-960","O-1020","Fe-1378",
"Fe-1438","Fe-1468","Fe-1229","Fe-1259","Fe-1289", "Fe-1319","Fe-1349","Fe-1379","Fe-1409","Fe-1469",
"Fe-1499","Fe-1410","Fe-1500","Fe-1440","p-1260","p-1467","p-1228","p-1258","p-1288","p-1348","p-1408",
"p-1438","p-1229","p-1289","p-1319","p-1379","p-1439","He-1496","Fe-1827","Fe-1857","Fe-1887","Fe-1917","O-1527",
"p-1558","Fe-1618","O-1618","He-1648","Fe-1677","Fe-1737"]
jobs21 = ["Fe-29","Fe-30","Fe-59","Fe-60","Fe-149","Fe-118","Fe-90","Fe-120","O-60","Fe-90","p-90"]


def getCORSIKALists(basePath,energyDir):
  corsikaList=sorted(glob.glob(basePath+str(energyDir)+"/*.bz2"))
  return corsikaList

def getCorsikaFiles(basePath,energyList,n1,n2):
  # energyList=sorted([f.name for f in os.scandir(basePath) if f.is_dir()])
  print("energyDir",energyList)
  corsikaFiles=[]
  for ienergy in energyList:
    print("ienergy",ienergy)
    corsikaList=getCORSIKALists(basePath,ienergy)
    # print("corsikaList",corsikaList)
    # corsikaFiles.append(corsikaList[:11])
    # corsikaFiles += corsikaList[45:50]
    corsikaFiles += corsikaList[n1:n2]
    # corsikaFiles += corsikaList
  # print("corsika files",corsikaFiles)
  return corsikaFiles

n=0
def submitToCondorFile(corsikaFile,primary):
  makeSubFile(corsikaFile,primary)
  corsikaID = str(corsikaFile).split("/")[-1]
  corsikaID = str(corsikaID).split(".")[0]
  # simulatedFile = simFiles+"dataSetUnique/"+primary+corsikaID+"GenDetFiltProcUnique.i3.gz"
  simulatedFile = simFiles+"dataSetGen1_6/"+primary+corsikaID+"Gen.i3.bz2"
  # simulatedFile = simFiles+"dataSetUnique/"+primary+corsikaID+"GenDetFiltProcUnique.i3.gz"
  print("simulated file",simulatedFile)
  if not os.path.exists(simulatedFile):
    corsikaID = int(''.join(i for i in corsikaID if i.isdigit()))
    batName = "{0}-{1}".format(primary,corsikaID)
    global n
    if batName not in queueJobs:
      subprocess.call(["condor_submit tempSubmit.sub -batch-name {0}-{1}".format(primary,corsikaID)], shell=True)
      print("batch name",batName)
      # if batName not in jobs14:
      if batName in jobs14:
        print("longer")
        # subprocess.call(["condor_submit tempSubmit.sub -batch-name {0}-{1}".format(primary,corsikaID)], shell=True)
      print('primary this has to be submitted lets see',primary,corsikaID,++n)
  subprocess.call(["rm tempSubmit.sub"], shell=True)
def submitToCondor(corsikaFiles,primary): 
  for ifiles in corsikaFiles:
    submitToCondorFile(ifiles,primary)
istep = 1
# for n1 in np.linspace(0,49,50):
# for n1 in np.linspace(50,99,50):
# for n1 in np.linspace(100,149,50):
for n1 in np.linspace(150,199,50):
  # # print("submitting proton files: ", basePathProton)
  # print("submitting Iron files from: ", basePathIron)
  ironCorsikaFiles = getCorsikaFiles(basePathIron,energyList,int(n1),int(n1+istep))
  submitToCondor(ironCorsikaFiles,"Fe")

  protonCorsikaFiles = getCorsikaFiles(basePathProton,energyList,int(n1),int(n1+istep))
  submitToCondor(protonCorsikaFiles,"p")

  heliumCorsikaFiles = getCorsikaFiles(basePathHelium,energyList,int(n1),int(n1+istep))
  submitToCondor(heliumCorsikaFiles,"He")

  oxygenCorsikaFiles = getCorsikaFiles(basePathOxygen,energyList,int(n1),int(n1+istep))
  submitToCondor(oxygenCorsikaFiles,"O")
  # #to add time in seconds +OriginalTime=21600