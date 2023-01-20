#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

PYTHON_SCRIPT=$BASEDIR/applyZenithRecoFilter.py
I3BUILD=/data/user/enpaudel/icecube_software/icetray_main/build
ICETRAY_ENV=$I3BUILD/env-shell.sh


$ICETRAY_ENV $PYTHON_SCRIPT 
echo "applyZenithRecoFilter.sh script was run successfully"