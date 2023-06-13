#!/bin/env python3

import os
import glob
import subprocess

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
	if energyID > 7.6:
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

# jobs7 = ["He-29","He-89","He-119","He-149","He-60","He-90","He-120","O-59","O-89","O-119","O-149","O-30","O-60","O-90","O-120","O-150"]
# jobs7 = ["p-30","p-60","p-90","p-120","p-150","p-88","p-89","He-29","He-89","He-119","He-149","He-60","He-90","He-120","He-28","He-58","He-88","He-117","He-56",
# "O-59","O-89","O-119","O-149","O-30","O-60","O-90","O-120","O-150","O-28","O-58","O-88","O-118","O-148","O-57","O-117",
# "Fe-29","Fe-30","Fe-59","Fe-60","Fe-89","Fe-90","Fe-119","Fe-120","Fe-149","Fe-150","Fe-58","Fe-118",
# "Fe-28","Fe-88","Fe-148","Fe-27","Fe-87","Fe-147","Fe-26"]
# jobs2 = ["p-145","p-175","p-625","p-715","p-775","p-805","p-865","p-895","He-505","He-535","He-565","He-625",
# "He-655","He-685","He-745","O-85","O-293","O-385","O-474","O-504","O-505","O-564","O-594","O-625","O-655","O-685"
# ,"O-715","O-745","O-714","O-2273","Fe-83","Fe-354","Fe-414","Fe-504","Fe-505","Fe-535","Fe-623","Fe-624",
# "Fe-625","Fe-654","Fe-744","Fe-803","Fe-833","Fe-1613","Fe-2153","Fe-2243","Fe-2783","Fe-2813"]
#energy 7.1 and 7.2
# jobs7 = ["p-3383","He-3443","He-3503","He-5962","O-263","O-443","O-593","O-923","O-504","O-594","O-1223","O-1313","O-1673","O-1973","O-2183","O-2573","O-2603",
# "O-2753","O-2963","O-3023","O-3053","O-3263","O-3443","O-3563","Fe-83","Fe-173","Fe-233","Fe-413","Fe-443","Fe-592","Fe-743","Fe-803","Fe-922","Fe-923",
# "Fe-1073","Fe-1103","Fe-1133","Fe-1732","Fe-1913","Fe-2063","Fe-2063","Fe-2093","Fe-2153","Fe-2243","Fe-2332","Fe-2423","Fe-2663","Fe-2693","Fe-2843","Fe-2902",
# "Fe-2903","Fe-3022","Fe-3083","Fe-3113","Fe-3172","Fe-3233","Fe-3263","Fe-3323","Fe-3352","Fe-3383","Fe-3442","Fe-3443","Fe-3502","Fe-3503","Fe-5933"]
# jobs7 = ["Fe-3773","Fe-4013","Fe-4103","Fe-4373","Fe-4613","Fe-4643","Fe-4762","Fe-5033","Fe-5273","Fe-5513"]
# jobs7 = ["Fe-6029","Fe-6030","Fe-6059","Fe-6060","p-6029","p-6030","p-6059","p-6060","O-6027","O-6028","O-6057","O-6058","He-6058","He-6059","He-6060","He-6028","He-6029","He-6030"]


#energy 7.3,7.4
# jobs7 = ["p-175","p-625","p-775","p-865","He-535","He-565","He-625","He-655","He-685","O-504","O-594","O-685","O-715","O-745","Fe-625","Fe-744"]
# jobs7 = ["Fe-3029","Fe-3030","He-3030","He-3060","O-3028","O-3058","Fe-3059","Fe-3060"]
# jobs7 = ["Fe-88","Fe-118","Fe-148","Fe-208","Fe-237","Fe-238","Fe-268","Fe-298","Fe-326","Fe-327","Fe-328","Fe-387","Fe-417","Fe-506","Fe-536","Fe-596","Fe-597",
# "p-238","p-420","p-450","p-508","p-538","p-539","He-207","O-3028","O-3058"]
# jobs7 = ["He-59","He-60","Fe-59","Fe-60","O-59","O-60","p-59","p-60","He-89","He-90","Fe-89","Fe-90","O-89","O-90","p-89","p-90"]
# queueJobs = runJobs+idleLines
# queueJobs = runJobs
# jobs7 = ["He-357","O-117","O-357","O-387","O-507","O-597","O-596","O-535","Fe-326","Fe-356","Fe-446","Fe-506","Fe-658","Fe-897","Fe-688","Fe-718","Fe-778","Fe-808"]
# jobs7 = ["Fe-208"]
# jobs7 = ["p-387","O-3777","O-3807","He-30","He-447","p-1078","p-1168","p-1288","p-1079","p-1229"]
# jobs7 = ["p-569","p-1289","p-1319","p-1349","p-1379"]
# jobs7 = ["Fe-30"]
# jobs7 = ["He-928","He-1378"]
# jobs7 = ["O-28", "O-58", "O-118", "O-148", "Fe-236", "Fe-237", "Fe-267", "Fe-297", "Fe-208", "Fe-238", "Fe-268", "Fe-298", "He-178", "He-208", "He-238", "He-268", "He-298", "O-297", "O-238", "O-268", "O-298", "Fe-596", "Fe-776", "Fe-806", "Fe-866", "Fe-956", "Fe-1046", "Fe-1286", "Fe-1406", "Fe-1496", "Fe-387", "Fe-417", "Fe-567", "Fe-597", "Fe-627", "Fe-687", "Fe-747", "Fe-807", "Fe-837", "Fe-867", "Fe-927", "Fe-1047"]
# jobs7 = ["Fe-1107","p-358","p-448","p-778","O-298", "p-1168", "p-1288", "p-1408", "p-1498", "He-777", "He-897", "He-957", "He-628", "He-658", "He-808", "He-838", "He-868", "He-928", "He-1048", "He-1078", "He-1198", "He-1228", "He-1288", "He-1318", "He-1378", "He-1498",]
# jobs14 = ["Fe-29", "Fe-59", "Fe-119", "Fe-149", "Fe-30", "Fe-60", "Fe-90",]
# jobs14 = ["O-117","Fe-147","O-1227","O-1257","Fe-2127"]
jobs14 = ["Fe-1046","Fe-956","Fe-806","Fe-776","Fe-596","O-2906","Fe-4436","Fe-4916",]
# jobs14 = ["Fe-1"]
jobs21 = ["Fe-29","Fe-59","Fe-30","Fe-60","Fe-149","He-119","p-150","p-120","p-90","Fe-150","Fe-359","Fe-479","Fe-448","Fe-478","Fe-568","Fe-419","Fe-449"]

jobs7LargerMemory = ["He-898","He-838"] #6 GB
# queueJobs = ["Fe-30"]
queueJobs = []
jobs7 = ["Fe-390","Fe-480","Fe-1134","Fe-1164","Fe-1194","Fe-1075","Fe-1195","Fe-1106",
"Fe-1196","Fe-1077","Fe-1107","Fe-1137","Fe-1167","Fe-1197","Fe-1078","Fe-1108","Fe-1138",
"Fe-1168","Fe-1198","Fe-1079","Fe-1109","Fe-1139","Fe-1169","Fe-1199","Fe-1080","Fe-1110",
"Fe-1140","Fe-1170","Fe-1200","p-2936","p-2966","He-2875","He-2876","He-2936","He-2966",
"He-2996","O-2904","O-2934","O-2875","O-2905","O-2965","O-2876","O-2906","O-2936","Fe-3113",
"Fe-3024","Fe-3055","Fe-3085","Fe-3115","Fe-3145","Fe-3026","Fe-3056","Fe-3086","Fe-3146",
"p-3024","p-3115","p-3056","p-3116","p-3146","He-3025","He-3145","He-3026","He-3056",
"He-3086","He-3116","He-3146","Fe-3295","Fe-3176","Fe-3206","Fe-3236","Fe-3266","Fe-3296",
"p-3174","p-3294","p-3176","p-3266","p-3296","He-3294","He-3175","He-3205","He-3235","He-3265",
"He-3176","He-3236","O-3204","O-3234","O-3264","O-3175","O-3205","O-3265","O-3295","O-3176",
"O-3236","O-3266","O-3296","Fe-3323","Fe-3383","Fe-3354","Fe-3384","Fe-3414","Fe-3325","Fe-3355",
"Fe-3385","Fe-3415","Fe-3445","Fe-3326","Fe-3356","Fe-3416","Fe-3446",
"O-3384", "O-3414", "O-3385", "O-3326", "O-3356", "O-3386", "O-3416", "O-3446",
 "Fe-3503", "Fe-3474", "Fe-3564", "Fe-3594", "Fe-3565", "Fe-3595", "Fe-3506", "Fe-3536", "Fe-3596",
 "p-3535","p-3565","p-3506","He-3503"]
jobs14 = ["O-57","Fe-87","He-87","O-117","Fe-237","O-237","Fe-267","Fe-297","Fe-627","Fe-687",
"Fe-747","Fe-628","Fe-658","Fe-1047","p-1047","p-1018","p-1050","He-928"]


def getCORSIKALists(basePath,energyDir):
	corsikaList=sorted(glob.glob(basePath+str(energyDir)+"/*.bz2"))
	return corsikaList

def getCorsikaFiles(basePath,energyList):
	# energyList=sorted([f.name for f in os.scandir(basePath) if f.is_dir()])
	print("energyDir",energyList)
	corsikaFiles=[]
	for ienergy in energyList:
		print("ienergy",ienergy)
		corsikaList=getCORSIKALists(basePath,ienergy)
		# print("corsikaList",corsikaList)
		# corsikaFiles.append(corsikaList[:11])
		corsikaFiles += corsikaList[190:200]
		# corsikaFiles += corsikaList
	# print("corsika files",corsikaFiles)
	return corsikaFiles

def submitToCondorFile(corsikaFile,primary):
	makeSubFile(corsikaFile,primary)
	corsikaID = str(corsikaFile).split("/")[-1]
	corsikaID = str(corsikaID).split(".")[0]
	# simulatedFile = simFiles+"dataSetUnique/"+primary+corsikaID+"GenDetFiltProcUnique.i3.gz"
	simulatedFile = simFiles+"dataSetGen1_6/"+primary+corsikaID+"Gen.i3.bz2"
	# simulatedFile = simFiles+"dataSetUnique/"+primary+corsikaID+"GenDetFiltProcUnique.i3.gz"
	# print("simulated file",simulatedFile)
	if not os.path.exists(simulatedFile):
		corsikaID = int(''.join(i for i in corsikaID if i.isdigit()))
		batName = "{0}-{1}".format(primary,corsikaID)
		if batName not in queueJobs:
			subprocess.call(["condor_submit tempSubmit.sub -batch-name {0}-{1}".format(primary,corsikaID)], shell=True)
			print("batch name",batName)
			# if batName not in jobs7:
			if batName in jobs14:
				print("longer")
				# subprocess.call(["condor_submit tempSubmit.sub -batch-name {0}-{1}".format(primary,corsikaID)], shell=True)
			print('primary',primary,corsikaID)
	subprocess.call(["rm tempSubmit.sub"], shell=True)
def submitToCondor(corsikaFiles,primary):	
	for ifiles in corsikaFiles:
		submitToCondorFile(ifiles,primary)

# # print("submitting proton files: ", basePathProton)
# print("submitting Iron files from: ", basePathIron)
ironCorsikaFiles = getCorsikaFiles(basePathIron,energyList)
submitToCondor(ironCorsikaFiles,"Fe")

protonCorsikaFiles = getCorsikaFiles(basePathProton,energyList)
submitToCondor(protonCorsikaFiles,"p")

heliumCorsikaFiles = getCorsikaFiles(basePathHelium,energyList)
submitToCondor(heliumCorsikaFiles,"He")

oxygenCorsikaFiles = getCorsikaFiles(basePathOxygen,energyList)
submitToCondor(oxygenCorsikaFiles,"O")
#to add time in seconds +OriginalTime=21600