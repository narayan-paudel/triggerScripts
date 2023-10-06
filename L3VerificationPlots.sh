#!/bin/bash
HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/setup.sh`

I3BASE=/data/user/enpaudel/icecube_software/icetray_work
# I3BASE=/data/user/enpaudel/icecube_software/icetray_current
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh

PYTHON_SCRIPT=$BASEDIR/L3VerificationPlots.py

fileList="$@"

$ICETRAY_ENV $PYTHON_SCRIPT -i $fileList --writeh5
echo "L3VerificationPlots.sh script was run successfully"