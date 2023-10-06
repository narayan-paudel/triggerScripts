#!/bin/bash

HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/


# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/setup.sh`
# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

# I3BASE=/data/user/enpaudel/icecube_software/icetray_main
I3BASE=/data/user/enpaudel/icecube_software/icetray_current
# I3BASE=/data/user/enpaudel/icecube_software/surfaceArray
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh


# GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/background/i3_daq/SPSSmokeTest/PFGCD_Run00138014_Subrun00000000.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HG.i3.gz
GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HGDomsets.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz
SEED=111

OUTPUT_FOLDER=/home/enpaudel/icecube/triggerStudy/simFilesTest/


PRIMARY_NAME=Fe
CORSIKA_ID=DAT002095

# PRIMARY_NAME=p
# CORSIKA_ID=DAT001164

# PRIMARY_NAME=Fe
# CORSIKA_ID=DAT000559


DETECTOR_PY=$I3SRC/simprod-scripts/resources/scripts/detector.py
FLAGS2="--UseGSLRNG --gcdfile ${GCD} --noInIce --LowMem --seed ${SEED} --nproc 1 --DetectorName IC86.2019 --no-FilterTrigger"
# FLAGS2+=" --inputfile $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2"
FLAGS2+=" --inputfile $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2"
FLAGS2+=" --output $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2"

echo running $DETECTOR_PY
echo flags $FLAGS2

$ICETRAY_ENV ${DETECTOR_PY} $FLAGS2

