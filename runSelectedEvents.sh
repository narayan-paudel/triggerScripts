#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

PYTHON_SCRIPT=$BASEDIR/selectedEvents.py
ICETRAY_ENV=/home/enpaudel/icecube/surfaceArray/build/env-shell.sh
# libraryAddress="/data/user/enpaudel/polarizationStudy/data/"

fileDir="/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/"
# fileList=$(ls -d $fileDir/DAT*GenDetFiltProc.i3.bz2)
fileList="$fileDir/DAT002192GenDetFiltProc.i3.bz2"
fileList="$fileList $fileDir/DAT002403GenDetFiltProc.i3.bz2"
fileList="$fileList $fileDir/DAT002289GenDetFiltProc.i3.bz2"
fileList="$fileList $fileDir/DAT002463GenDetFiltProc.i3.bz2"
fileList="$fileList $fileDir/DAT002218GenDetFiltProc.i3.bz2"
# fileList=${$fileList[@]:3:10}
echo "file list" $fileList
arr=($fileList)
echo "starting with first shower" ${arr[0]}
$ICETRAY_ENV $PYTHON_SCRIPT $fileList
# $ICETRAY_ENV $PYTHON_SCRIPT "${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]} ${arr[4]}"
echo "runTrigcheck script was run successfully"