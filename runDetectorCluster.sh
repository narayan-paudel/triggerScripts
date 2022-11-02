#!/bin/bash

HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/


eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

I3BASE=/data/user/enpaudel/icecube_software/icetray_main
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh

# ICETRAY_ENV=/home/enpaudel/icecube/icetray/build/env-shell.sh
# I3BUILD=/home/enpaudel/icecube/icetray/build/
# I3SRC=/home/enpaudel/icecube/icetray/src/

# GCD=/data/user/kath/testdata/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305Updated.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/GCD-Survey-AntITScint_2020.02.24_Snow210305.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GCD-Survey-AntITScint_2020.02.24_newCalibration.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GCD-Survey-AntITScint_2020.02.24_newDetectorStatus.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305OldDetectorStatus.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz
# GCD="/home/enpaudel/icecube/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305ModifiedNoSMTDOMSet.i3.gz"

######################################################################################
# GCD=/home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz
# GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz
GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz
######################################################################################


# $ICETRAY_ENV $PYTHON_SCRIPT 
# /home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz
#input ex:/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/iron/lgE_17.0/sin2_0.9/000054/DAT000054
#input ex:/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.2/DAT001203.bz2"
#input ex:/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.7/DAT000638.bz2"
#input ex:/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/6.0/DAT000671.bz2"

SCRATCH_FOLDER=/scratch/enpaudel/
CORSIKA_FILE_PATH=$1
CORSIKA_FOLDER=$(dirname $(realpath -s $1))
CORSIKA_FOLDER=$CORSIKA_FOLDER/
CORSIKA_ID_BZ=$(basename $1)
CORSIKA_ID=${CORSIKA_ID_BZ%.*}
corsikaRunID=$(echo $CORSIKA_ID | sed 's/[^0-9]*//g'|sed 's/^0*//')
echo using runID $corsikaRunID
logFile=$CORSIKA_FOLDER$CORSIKA_ID.log
if [[ -s "$logFile" ]]; then
	ERANGE=$(sed -n -e '/PRIMARY ENERGY =/p' $logFile)
	ENERGY=$(echo $ERANGE | sed -e 's/[^0-9.]//g')
	PRIMARY=$(sed -n -e '/PRMPAR/p' $logFile)
	PRIMARY=$(echo $PRIMARY | sed -e 's/[^0-9]//g')
	if [[ $PRIMARY == 14 ]]; then
		PRIMARY_NAME=p
		pSEED=14
	elif [[ $PRIMARY == 5626 ]]; then
		PRIMARY_NAME=Fe
		pSEED=5626
	elif [[ $PRIMARY == 402 ]]; then
		PRIMARY_NAME=He
		pSEED=402
	elif [[ $PRIMARY == 1608 ]]; then
		PRIMARY_NAME=O
		pSEED=1608
	fi
	echo log file exists
	discR=$(python -c "import math;print(800+(int(math.log10($ENERGY))-5)*300 + (2*(int(math.log10($ENERGY))-5)//3)*300 + ((int(math.log10($ENERGY))-5)//3)*300)")
else
	PRIMARY=$(basename $(dirname $CORSIKA_FOLDER))
	ENERGY=$(basename $CORSIKA_FOLDER)
	if [[ $PRIMARY == 10410 ]]; then
		PRIMARY_NAME=p
		pSEED=14
	elif [[ $PRIMARY == 10889 ]]; then
		PRIMARY_NAME=Fe
		pSEED=5626
	elif [[ $PRIMARY == 11663 ]]; then
		PRIMARY_NAME=He
		pSEED=402
	elif [[ $PRIMARY == 12605 ]]; then
		PRIMARY_NAME=O
		pSEED=1608
	fi
	discR=$(python -c "import math;print(800+(int($ENERGY)-5)*300 + (2*(int($ENERGY)-5)//3)*300 + ((int($ENERGY)-5)//3)*300)")
fi
echo using disc radius $discR
echo energy $ENERGY
nSamples=100
echo nSamples $nSamples

OUTPUT_FOLDER=/home/enpaudel/icecube/triggerStudy/simFiles/
# OUTPUT_FOLDER=/home/enpaudel/icecube/triggerStudy/simFilesTest/
start_time=$SECONDS
# echo start time $start_time
###########################################################################
# bunzip2 -ckd $1 > $OUTPUT_FOLDER/temp/${PRIMARY_NAME}$CORSIKA_ID

if [ ! -f $OUTPUT_FOLDER/temp/${PRIMARY_NAME}$CORSIKA_ID ]; then
	bunzip2 -ckd $1 > $OUTPUT_FOLDER/temp/${PRIMARY_NAME}$CORSIKA_ID
fi
# time_zip=$(($SECONDS - $start_time))
##########################################################################
#calculate the seed
echo calculating the seed
echo $corsikaRunID
echo $PRIMARY_NAME
echo $pSEED
echo $ENERGY
eSEED=$(python -c "import math;print(int(math.log10($ENERGY)))")
SEED=$(($corsikaRunID+($pSEED*10+$eSEED)*10800000))
# eSEED=$(python -c "import math;print(int(math.log10($ENERGY)*10))")
# SEED=$(($corsikaRunID+($pSEED*10+$eSEED)*108000000)) #
echo $SEED

GENERATOR_PY=$I3SRC/simprod-scripts/resources/scripts/icetopshowergenerator.py 

FLAGS="--UseGSLRNG --gcdfile ${GCD} --seed ${SEED} --nproc 1 --RunID $corsikaRunID --samples $nSamples --r $discR --no-PropagateMuons"
FLAGS="$FLAGS --inputfile $OUTPUT_FOLDER/temp/${PRIMARY_NAME}$CORSIKA_ID"
FLAGS="$FLAGS --output $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenTemp.i3.bz2"

echo running $GENERATOR_PY
echo flags $FLAGS
# $ICETRAY_ENV ${GENERATOR_PY} $FLAGS
####################################################
$ICETRAY_ENV ${GENERATOR_PY} $FLAGS

# if [ ! -f $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2 ]; then
# 	$ICETRAY_ENV ${GENERATOR_PY} $FLAGS
# fi

time_gen=$(($SECONDS - $start_time))
##########################################################################
rm $OUTPUT_FOLDER/temp/${PRIMARY_NAME}$CORSIKA_ID
mv $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenTemp.i3.bz2 $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2

# DETECTOR_PY=$I3SRC/simprod-scripts/resources/scripts/detector.py
# FLAGS2="--UseGSLRNG --gcdfile ${GCD} --noInIce --LowMem --seed ${SEED} --nproc 1 --DetectorName IC86.2019 --no-FilterTrigger"
# FLAGS2+=" --inputfile $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2"
# FLAGS2+=" --output $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2"

# echo running $DETECTOR_PY
# echo flags $FLAGS2

# $ICETRAY_ENV ${DETECTOR_PY} $FLAGS2

# PHOTONDIR="/cvmfs/icecube.opensciencegrid.org/data/photon-tables"
# FILTERING_PY=$I3SRC/filterscripts/resources/scripts/SimulationFiltering.py
# FLAGS3=" --needs_wavedeform_spe_corr --photonicsdir ${PHOTONDIR} -g ${GCD}"
# FLAGS3+=" --input $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2"
# FLAGS3+=" --output $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenDetFilt.i3.bz2"

# echo running $FILTERING_PY
# echo flags $FLAGS3

# $ICETRAY_ENV ${FILTERING_PY} $FLAGS3
# # if [ ! -f $OUTPUT_FOLDER/dataSet/${CORSIKA_ID}GenDetFilt.i3.bz2 ]; then
# # 	$ICETRAY_ENV ${FILTERING_PY} $FLAGS3
# # fi

# # time_filt=$(($SECONDS - $start_time))

# PROCESS_PY=$I3SRC/filterscripts/resources/scripts/offlineL2/process.py
# # PROCESS_PY=$I3SRC/filterscripts/python/pass2/processing.py
# FLAGS4=" -s --photonicsdir ${PHOTONDIR} -g ${GCD}"
# FLAGS4+=" --input $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenDetFilt.i3.bz2"
# FLAGS4+=" --output $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}GenDetFiltProc.i3.bz2"

# echo running $PROCESS_PY
# echo flags $FLAGS4

# echo "$ICETRAY_ENV ${PROCESS_PY} $FLAGS4"
# $ICETRAY_ENV ${PROCESS_PY} $FLAGS4





# ###########################################################################################################
# ###########################################################################################################
# # GENERATOR_PY=$I3BUILD/simprod-scripts/resources/scripts/icetopshowergenerator.py --UseGSLRNG\
# #  --gcdfile /home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz --seed 0 --nproc 1 --samples 100 --inputfile\
# #   /home/enpaudel/icecube/triggerStudy/scripts/DAT059871 --output /home/enpaudel/icecube/triggerStudy/detectorResponse/test.i3.bz

# # DETECTOR_PY=$I3build/simprod-scripts/resources/scripts/detector.py --output /home/enpaudel/icecube/triggerStudy/detectorResponse/test2019detectorNoFilter.i3.bz\
# #   --inputfile /home/enpaudel/icecube/triggerStudy/detectorResponse/test2019.i3.bz --UseGSLRNG --gcdfile /home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz\
# #    --noInIce --LowMem --seed 0 --nproc 0 --DetectorName IC86.2019 --no-FilterTrigger

# # GENERATOR_PY=$I3BUILD/simprod-scripts/resources/scripts/icetopshowergenerator.py\
# # --UseGSLRNG --gcdfile /data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz\
# # --seed 0 --nproc 0 --samples 100 --inputfile data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.0/DAT000001.bz2\
# # --output /home/enpaudel/icecube/triggerStudy/detectorResponse/test.i3.bz
# # DETECTOR_PY=$I3BUILD/simprod-scripts/resources/scripts/detector.py --inputfile data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.0/DAT000001.bz2\
# #  --output /home/enpaudel/icecube/triggerStudy/detectorResponse/test.i3.bz --UseGSLRNG --gcdfile /data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz\
# #   --noInIce --LowMem --seed 0 --nproc 0 --DetectorName IC86.2012

# #To use directly on command line:
# # /home/enpaudel/icecube/surfaceArray/build/simprod-scripts/resources/scripts/icetopshowergenerator.py --UseGSLRNG --gcdfile /home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz --seed 0 --nproc 1 --samples 100 --inputfile /home/enpaudel/icecube/triggerStudy//simFiles/temp/DAT059871 --output /home/enpaudel/icecube/triggerStudy/detectorResponse/test2019.i3.bz
# # /home/enpaudel/icecube/surfaceArray/build/simprod-scripts/resources/scripts/detector.py --output /home/enpaudel/icecube/triggerStudy/detectorResponse/test2019detector.i3.bz  --inputfile /home/enpaudel/icecube/triggerStudy/detectorResponse/test2019.i3.bz --UseGSLRNG --gcdfile /home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz --noInIce --LowMem --seed 0 --nproc 0 --DetectorName IC86.2019 --no-FilterTrigger
# # /home/enpaudel/icecube/surfaceArray/build/simprod-scripts/resources/scripts/detector.py --output /home/enpaudel/icecube/triggerStudy/detectorResponse/test2019detectorNoFilter.i3.bz  --inputfile /home/enpaudel/icecube/triggerStudy/detectorResponse/test2019.i3.bz --UseGSLRNG --gcdfile /home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz --noInIce --LowMem --seed 0 --nproc 0 --DetectorName IC86.2019
