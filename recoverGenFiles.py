#!/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"

pathUnique = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique/"
pathGen = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetGen/"

UniqueList = os.listdir(pathUnique)
UniqueList = [f[:-23] for f in UniqueList]
GenList = os.listdir(pathGen)
GenList = [f[:-7] for f in GenList]
UniqueGenList = [f for f in GenList if f in UniqueList]
# UniqueGenList = [f+".i3.bz2" for f in UniqueGenList if f not in UniqueList]
# extraUniqueList = [f for f in UniqueList if f not in GenList]
# print(extraUniqueList)
print(len(GenList))
print(len(UniqueList))
print(len(UniqueGenList))