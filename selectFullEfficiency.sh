#!/bin/bash

HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/


eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

# I3BASE=/data/user/enpaudel/icecube_software/icetray_main
I3BASE=/data/user/enpaudel/icecube_software/surfaceArray
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh


SELECT_PY=/home/enpaudel/icecube/triggerStudy/triggerScripts/selectFullEfficiency.py
INPUT=$1



$ICETRAY_ENV ${SELECT_PY} $INPUT
#############################################
#############################################