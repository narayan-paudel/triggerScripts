#!/bin/env python3

import numpy as np

import glob
import re


def energyToBin(E):
  return int(round((E-5.0)*10))


def corsikaListEnergy(lowE,highE):
  lcol = energyToBin(lowE)
  hcol = energyToBin(highE)
  longList = np.linspace(1,5971,200)
  longList = [int(n) for n in longList]
  newList = []
  for i in range(lcol,hcol+1):
    incrementList = [x+i for x in longList]
    # print(incrementList)
    newList += incrementList
  newList=sorted(newList)
  return newList

corsikaList = corsikaListEnergy(7.0,7.5)
# print(corsikaList)

listOFCorsikaFiles = []

# uniqueFiles = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique/pDAT005995GenDetFiltProcUnique.i3.gz"
uniqueFilesPath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique/"



uniqueFiles = sorted(glob.glob(uniqueFilesPath+"/*"))
for ifile in uniqueFiles:
  corsikaID = str(ifile).split("/")[-1]
  corsikaID = str(corsikaID).split(".")[0]
  corsikaID = int(re.findall(r'\d+',corsikaID)[0])
  if corsikaID in corsikaList:
    print(ifile)



