import numpy as np
import re
import subprocess

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
idleLines = idleJobs("../idleJobs.txt")

def delJobs(idleJobs):
	for ijobs in idleJobs:
		subprocess.call(["condor_rm {}".format(ijobs)], shell=True)

delJobs(idleLines)
