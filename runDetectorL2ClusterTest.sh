#!/bin/bash

HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/


eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

I3BASE=/data/user/enpaudel/icecube_software/icetray_main
# I3BASE=/data/user/enpaudel/icecube_software/surfaceArray
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh


GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz
SEED=101

OUTPUT_FOLDER=/home/enpaudel/icecube/triggerStudy/simFilesTest/


# PRIMARY_NAME=Fe
# CORSIKA_ID=DAT002095

PRIMARY_NAME=p
CORSIKA_ID=DAT001164

DETECTOR_PY=$I3SRC/simprod-scripts/resources/scripts/detector.py
FLAGS2="--UseGSLRNG --gcdfile ${GCD} --noInIce --LowMem --seed ${SEED} --nproc 1 --DetectorName IC86.2019 --no-FilterTrigger"
# FLAGS2+=" --inputfile $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2"
FLAGS2+=" --inputfile $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2"
FLAGS2+=" --output $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2"

echo running $DETECTOR_PY
echo flags $FLAGS2

$ICETRAY_ENV ${DETECTOR_PY} $FLAGS2
# if [ ! -f $OUTPUT_FOLDER/dataSet/${CORSIKA_ID}GenDet.i3.bz2 ]; then
# 	$ICETRAY_ENV ${DETECTOR_PY} $FLAGS2
# fi
# time_det=$(($SECONDS - $start_time))


PHOTONDIR="/cvmfs/icecube.opensciencegrid.org/data/photon-tables"
FILTERING_PY=$I3SRC/filterscripts/resources/scripts/SimulationFiltering.py
FLAGS3=" --needs_wavedeform_spe_corr --photonicsdir ${PHOTONDIR} -g ${GCD}"
FLAGS3+=" --input $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2"
FLAGS3+=" --output $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}GenDetFilt.i3.bz2"

echo running $FILTERING_PY
echo flags $FLAGS3

$ICETRAY_ENV ${FILTERING_PY} $FLAGS3
# # if [ ! -f $OUTPUT_FOLDER/dataSet/${CORSIKA_ID}GenDetFilt.i3.bz2 ]; then
# # 	$ICETRAY_ENV ${FILTERING_PY} $FLAGS3
# # fi

# time_filt=$(($SECONDS - $start_time))
#####################################################use either process.py(baseprocessing.py) or processing.py####
PROCESS_PY=$I3SRC/filterscripts/resources/scripts/offlineL2/process.py
FLAGS4=" -s --photonicsdir ${PHOTONDIR} -g ${GCD}"
FLAGS4+=" --input $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}GenDetFilt.i3.bz2"
FLAGS4+=" --output $OUTPUT_FOLDER/dataSetGen/${PRIMARY_NAME}${CORSIKA_ID}GenDetFiltProc.i3.bz2"

echo running $PROCESS_PY
echo flags $FLAGS4

echo "$ICETRAY_ENV ${PROCESS_PY} $FLAGS4"
$ICETRAY_ENV ${PROCESS_PY} $FLAGS4
############################################################################

# # echo "time_diff=[$SECONDS - $start_time]"
# time_proc=$(($SECONDS - $start_time))
# # echo $time_diff
# # printf "$ENERGY $time_zip $time_gen $time_det $time_filt $time_proc \n" >> $BASEDIR/../energyTimeMulti.txt  
# printf "$ENERGY $time_proc \n" >> $BASEDIR/../energyTimeMultiStudyProper.txt
# echo removing intermediate files
# ################################# rm $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2
# rm $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenDetFilt.i3.bz2 $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2 

# UNIQUE_PY=$BASEDIR/addUniqueRunID.py
# INPUT_FILE=$OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenDetFiltProc.i3.bz2

# $ICETRAY_ENV $UNIQUE_PY $INPUT_FILE
# # rm $ifile
