import numpy as np
import re
import subprocess

# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('lcol', type=int, default=0, help="starting columns for counting")
# parser.add_argument('hcol', type=int, default=30, help="upper columns for counting")
# args = parser.parse_args()

def readTxt(filePath):
	with open(filePath,"r") as f:
		lines = f.readlines()
	lines = [iline.split("/")[-1] for iline in lines ]
	lines = [iline.split(".")[0] for iline in lines ]
	regex = re.compile(r"\d+")
	lines = [int(x) for iline in lines for x in regex.findall(iline)]
	print("lines",lines)

# readTxt("./processedHe.txt")


def runningJobs(filePath):
	with open(filePath,"r") as f:
		next(f)
		next(f)
		lines = f.readlines()
	lines = [iline.split(" ")[1] for iline in lines if iline.split(" ")[0] == "enpaudel" ]
	print("lines",lines)
	# lines = [iline.split("/")[-1] for iline in lines ]
	# lines = [iline.split(".")[0] for iline in lines ]
	# regex = re.compile(r"\d+")
	# lines = [int(x) for iline in lines for x in regex.findall(iline)]
	# print("lines",lines)
	return lines
# lines = runningJobs("../runningJobs.txt")


def idleJobs(filePath):
	with open(filePath,"r") as f:
		next(f)
		next(f)
		lines = f.readlines()
	print("Hello hello")
	print("lines",lines[:3])
	lines = [iline.split(" ")[-1].rstrip() for iline in lines if iline.split(" ")[0] == "enpaudel" ]
	# lines = [iline.split("/\")[-1] for iline in lines if iline.split(" ")[0] == "enpaudel" ]
	print("lines",lines[:3])
	# lines = [iline.split("/")[-1] for iline in lines ]
	# lines = [iline.split(".")[0] for iline in lines ]
	# regex = re.compile(r"\d+")
	# lines = [int(x) for iline in lines for x in regex.findall(iline)]
	# print("lines",lines)
	return lines
# idleLines = idleJobs("../idleJobs.txt")
# print("idle lines",idleLines)
def doubleJobs(filePath):
	with open(filePath,"r") as f:
		next(f)
		next(f)
		next(f)
		next(f)
		lines = f.readlines()
	print("Hello hello")
	lines = lines[:-5]
	print("lines",lines[0])
	# print([line.split(" ")[-4] for line in lines[0:238]])
	# lines = [iline.split(" ")[-1].rstrip() for iline in lines if iline.split(" ")[1] == "enpaudel" ]
	lines = [iline.split(" ")[-1].rstrip() for iline in lines if iline.split(" ")[-4] == "2" ]
	# lines1 = [lines[n].split(" ")[-4] for n,iline in enumerate(lines)]
	# lines = [iline.split(" ")[1][0] for iline in lines]
	# lines = [iline.split("/\")[-1] for iline in lines if iline.split(" ")[0] == "enpaudel" ]
	# print("lines",lines)
	# lines = [iline.split("/")[-1] for iline in lines ]
	# lines = [iline.split(".")[0] for iline in lines ]
	# regex = re.compile(r"\d+")
	# lines = [int(x) for iline in lines for x in regex.findall(iline)]
	# print("lines",lines)
	lines = [float(iline) for iline in lines]
	return lines
idleLines = doubleJobs("../idleJobs.txt")
print("double lines",idleLines, len(idleLines))

def delJobs(idleJobs):
	for ijobs in idleJobs:
		subprocess.call(["condor_rm {}".format(ijobs)], shell=True)

###################### delJobs(idleLines)
#need to run cmd ls -l  ../simFiles/dataSet/He*Proc* > ../completedJobs.txt
# ls -l  ../simFiles/dataSetUnique1_6/p*Proc* > ../completedJobs.txt
################################################3
# counting completed jobs
################################################

# print(longList)

###############################################
##############################################
# longList = np.linspace(1,5971,200)
# # longList = np.linspace(1,3001,101)
# longList = [int(n) for n in longList]
# # print("longList",longList)
# newList = []
# for i in range(23,24):
# 	incrementList = [x+i for x in longList]
# 	# print(incrementList)
# 	newList += incrementList
# newList=sorted(newList)
# print("longList",newList)

def completedJobs(filePath,newList):
	with open(filePath,"r") as f:
		# next(f)
		# next(f)
		lines = f.readlines()
	# print("Hello hello")
	lines = [iline.split(" ")[-1].rstrip() for iline in lines if iline.split(" ")[2] == "enpaudel" ]
	lines = [iline.split("/")[-1] for iline in lines ]
	lines = [iline.split(".")[0] for iline in lines ]
	regex = re.compile(r"\d+")
	lines = [int(x) for iline in lines for x in regex.findall(iline)]
	# print("lines",lines)
	remainElt = [x for x in newList if x not in lines]
	# print("remainElt",remainElt)
	# print("lines",lines)

	return lines


def checkCompletedFiles(lcol,hcol):
  longList = np.linspace(1,5971,200)
  longList = [int(n) for n in longList]
  newList = []
  for i in range(lcol,hcol):
    incrementList = [x+i for x in longList]
    # print(incrementList)
    newList += incrementList
  newList=sorted(newList)
  cmpletedLines = completedJobs("../completedJobs.txt",newList)
  # print("completed jobs",len([i for i in newList if i in cmpletedLines]),[i for i in newList if i in cmpletedLines],len([i for i in newList if i in cmpletedLines]))
  print("completed jobs",hcol,len([i for i in newList if i in cmpletedLines]),len([i for i in newList if i in cmpletedLines]))



for i in range(0,30):
  checkCompletedFiles(i,i+1)

checkCompletedFiles(0,31)