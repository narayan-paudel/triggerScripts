#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

PYTHON_SCRIPT=$BASEDIR/addUniqueRunID.py
ICETRAY_ENV=/home/enpaudel/icecube/surfaceArray/build/env-shell.sh
# libraryAddress="/data/user/enpaudel/polarizationStudy/data/"

fileDir="/home/enpaudel/icecube/triggerStudy/simFiles/dataSet/"
# fileList=$(ls -d $fileDir/*DAT*GenDetFiltProc.i3.bz2)
fileList="$@"
echo filelist $fileList
for ifile in $fileList;do
	echo "starting with first shower" $ICETRAY_ENV $PYTHON_SCRIPT  $ifile
	$ICETRAY_ENV $PYTHON_SCRIPT $ifile
	echo "deleting the original file"
	rm $ifile
done
echo "unique run ID was added"