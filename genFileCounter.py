#!/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
basePathGen = "/data/user/enpaudel/triggerStudy/simFiles/dataSetGen/"
pGenFiles = sorted(glob.glob(basePathGen+"p*Gen.i3.bz2"))
# print("pgenfiles",pGenFiles)
for ifile in pGenFiles:
	subprocess.call(["python frameCounter.py {}".format(ifile)], shell=True)
