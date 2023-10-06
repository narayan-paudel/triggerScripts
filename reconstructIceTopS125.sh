#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/setup.sh`

I3BASE=/data/user/enpaudel/icecube_software/icetray_work
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh

PYTHON_SCRIPT=$BASEDIR/reconstructIceTopS125.py

fileList="$@"

$ICETRAY_ENV $PYTHON_SCRIPT -i $fileList
echo "runFilteredEvents.sh script was run successfully"