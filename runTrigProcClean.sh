#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

PYTHON_SCRIPT=$BASEDIR/triggerProcClean.py
I3BASE=/data/user/enpaudel/icecube_software/surfaceArray
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh
# libraryAddress="/data/user/enpaudel/polarizationStudy/data/"

fileDir="/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/"
fileList=$(ls -d $fileDir/*DAT*GenDetFiltProc.i3.bz2)
# fileList=${$fileList[@]:3:10}
arr=($fileList)
echo "starting with first shower" $ICETRAY_ENV $PYTHON_SCRIPT  ${arr[0]}
$ICETRAY_ENV $PYTHON_SCRIPT $fileList
echo "runTrigcheck script was run successfully"