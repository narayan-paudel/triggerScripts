#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

PYTHON_SCRIPT=$BASEDIR/triggerCheck.py
ICETRAY_ENV=/home/enpaudel/icecube/surfaceArray/build/env-shell.sh
# libraryAddress="/data/user/enpaudel/polarizationStudy/data/"

fileDir="/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/"
fileList=$(ls -d $fileDir/DAT*GenDetFiltProc.i3.bz2)
for ifile in $fileList;do
	echo "starting with first shower" $ICETRAY_ENV $PYTHON_SCRIPT  $ifile
	$ICETRAY_ENV $PYTHON_SCRIPT $ifile
done
echo "runTrigcheck script was run successfully"