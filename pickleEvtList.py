#!/usr/bin/env python3


from customColors import qualitative_colors

import numpy as np
import pickle


from inclinedTriggerTools import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, nargs='+', default="/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.hdf5", help='Input files after running detector.py.')
args = parser.parse_args()


basePath = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/"

plotFolder = "/home/enpaudel/icecube/triggerStudy/plots/"

#example python pickleEvtList.py /home/enpaudel/icecube/triggerStudy/simFiles/dataSetCleanTanks/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.hdf5


triggerList = ["HLC6_5000","tank6_5000","tank6_4000","tank6_3000","tank6_2000",
  "tank7_5000","tank7_4000","tank7_3000","tank7_2000","tank8_5000","tank8_4000","tank8_3000","tank8_2000"]
triggerList7 = ["HLC6_5000","tank7_5000","tank7_4000","tank7_3000","tank7_2000"]
triggerListSelect = ["HLC6_5000","tank6_3000","tank7_3000","tank8_3000"]


trigWindow = 10**(-6) # in ns

sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
# sin2ZenBins = [0.0,0.822]

energyBins = 10**np.linspace(5, 8.0, 31)
energyBinsShort = 10**np.linspace(5, 8.0, 7)
# energyBins = 10**np.linspace(5, 8.0, 7)
energyBinslgE = np.linspace(5.0,8.9,4000)
energyBinCenter = [5.1,6.1,7.1,8.1]

print("args.input",args.input)

pickleFile = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetPickle1_6/{}.pkl".format(args.input[0].split("/")[-1].split(".")[0])

print("picklefile",pickleFile)

# recoI3 = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetInclFilt/showerInclinedFilter.hdf5"

  
evtList = extractEvents(args.input)
with open(pickleFile,"wb") as f:
  pickle.dump(evtList,f)


