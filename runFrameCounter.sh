#!/bin/bash

# fileList=/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/*Filt*
fileList=/home/enpaudel/icecube/triggerStudy/simFiles/dataSet/*Gen.i3.bz2
# fileList=/home/enpaudel/icecube/triggerStudy/simFiles/dataSetTest/testFile/*Gen.i3.bz2

for ifile in $fileList;do
	# echo "counting" $ifile
	python frameCounter.py $ifile
done
