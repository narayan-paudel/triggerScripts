#!/bin/bash

# fileList=/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/*Filt*
fileList=/home/enpaudel/icecube/triggerStudy/simFiles/trigChk/*i3*

for ifile in $fileList;do
	echo "counting" $ifile
	python frameCounter.py $ifile
done
