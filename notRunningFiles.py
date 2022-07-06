import numpy as np
import re
import subprocess

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('lcol', type=int, default=0, help="starting columns for counting")
parser.add_argument('hcol', type=int, default=30, help="upper columns for counting")
args = parser.parse_args()

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
	# print("lines",lines)
	# lines = [iline.split("/")[-1] for iline in lines ]
	# lines = [iline.split(".")[0] for iline in lines ]
	# regex = re.compile(r"\d+")
	# lines = [int(x) for iline in lines for x in regex.findall(iline)]
	# print("lines",lines)
	return lines
lines = runningJobs("../runningJobs.txt")
# print("sorted",sorted(lines,key=lambda x:int(re.search(r'\d+', x).group())))
def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
for x in sorted_nicely(lines):
    print(x)



################################################
longList = np.linspace(1,5971,200)
longList = [int(n) for n in longList]
# print(longList)
newList = []
for i in range(args.lcol,args.hcol):
	incrementList = [x+i for x in longList]
	# print(incrementList)
	newList += incrementList
newList=sorted(newList)




remainElt = [x for x in newList if x not in lines]
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


# print("not running jobs",remainElt)
