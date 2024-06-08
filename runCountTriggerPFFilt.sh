#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

############################
# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
# I3BASE=/data/user/enpaudel/icecube_software/icetray_current
############################

# I3BASE=/data/user/enpaudel/icecube_software/icetray_main
# I3BASE=/data/user/enpaudel/icecube_software/surfaceArray
# I3BASE=/data/user/enpaudel/icecube_software/icetray_current
# I3BASE=/data/user/enpaudel/icecube_software/icetray_work

####################################
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/setup.sh`
I3BASE=/data/user/enpaudel/icecube_software/icetray_current
####################################
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh

# PYTHON_SCRIPT=$BASEDIR/countTriggerPFFilt.py
PYTHON_SCRIPT=$BASEDIR/countTriggerEventPFFilt.py


# fileDir="/home/enpaudel/icecube/triggerStudy/simFiles/ITGen/"
# fileList=$(ls -d $fileDir/??DAT*GenDetFiltProc.i3.bz2)

$ICETRAY_ENV $PYTHON_SCRIPT -g $1 -i ${@:2}

#./runTrigProcCleanChunk.sh ../simFiles/dataSetUnique/FeDAT005775GenDetFiltProcUnique.i3.gz