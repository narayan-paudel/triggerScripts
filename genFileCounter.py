#!/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
# basePathGen = "/data/user/enpaudel/triggerStudy/simFiles/dataSetGen1_6/"
basePathGen = "/data/user/enpaudel/triggerStudy/simFiles/dataSet/"
pGenFiles = sorted(glob.glob(basePathGen+"*Gen.i3.bz2"))
# print("pgenfiles",pGenFiles)
for ifile in pGenFiles:
	subprocess.call(["python frameCounter.py {}".format(ifile)], shell=True)
