#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

PYTHON_SCRIPT=$BASEDIR/addUniqueRunID.py
I3BUILD=/data/user/enpaudel/icecube_software/icetray_main/build
ICETRAY_ENV=$I3BUILD/env-shell.sh
# libraryAddress="/data/user/enpaudel/polarizationStudy/data/"

fileDir="/home/enpaudel/icecube/triggerStudy/simFiles/dataSet/"
# fileDir="/home/enpaudel/icecube/triggerStudy/simFilesTest/dataSet/"
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