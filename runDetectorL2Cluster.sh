#!/bin/bash

HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/


eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

# I3BASE=/data/user/enpaudel/icecube_software/icetray_main
I3BASE=/data/user/enpaudel/icecube_software/icetray_inclinedTrigger
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
# GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz
######################################################################################


# $ICETRAY_ENV $PYTHON_SCRIPT 
# /home/acoleman/work/datasets/gcd-files/GCD-Survey-AntITScint_2020.02.24.i3.gz
#input ex:/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/iron/lgE_17.0/sin2_0.9/000054/DAT000054
#input ex:/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.2/DAT001203.bz2"
#input ex:/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.7/DAT000638.bz2"
#input ex:/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/6.0/DAT000671.bz2"
#input ex:/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.0/DAT000001.bz2"

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
	OBSLEV=$(sed -n -e '/OBSLEV/p' $logFile)
	OBSLEV=$(echo $OBSLEV | sed -e 's/[^0-9]//g')
	OBSLEV=${OBSLEV:0:4}
	if [[ $PRIMARY == 14 ]]; then
		PRIMARY_NAME=p
		# pSEED=14
		pSEED=1
	elif [[ $PRIMARY == 5626 ]]; then
		PRIMARY_NAME=Fe
		# pSEED=5626
		pSEED=4
	elif [[ $PRIMARY == 402 ]]; then
		PRIMARY_NAME=He
		# pSEED=402
		pSEED=2
	elif [[ $PRIMARY == 1608 ]]; then
		PRIMARY_NAME=O
		# pSEED=1608
		pSEED=3
	fi
	echo log file exists
	discR=$(python -c "import math;print(800+(int(math.log10($ENERGY))-5)*300 + (2*(int(math.log10($ENERGY))-5)//3)*300 + ((int(math.log10($ENERGY))-5)//3)*300)")
else
	PRIMARY=$(basename $(dirname $CORSIKA_FOLDER))
	ENERGY=$(basename $CORSIKA_FOLDER)
	if [[ $PRIMARY == 10410 ]]; then
		PRIMARY_NAME=p
		# pSEED=14
		pSEED=1
		OBSLEV=2834
	elif [[ $PRIMARY == 10889 ]]; then
		PRIMARY_NAME=Fe
		# pSEED=5626
		pSEED=4
		OBSLEV=2834
	elif [[ $PRIMARY == 11663 ]]; then
		PRIMARY_NAME=He
		# pSEED=402
		pSEED=2
		OBSLEV=2837
	elif [[ $PRIMARY == 12605 ]]; then
		PRIMARY_NAME=O
		# pSEED=1608
		pSEED=3
		OBSLEV=2837
	fi
	discR=$(python -c "import math;print(800+(int($ENERGY)-5)*300 + (2*(int($ENERGY)-5)//3)*300 + ((int($ENERGY)-5)//3)*300)")
fi

OBSLEV_GOOD=2840
RAISE_H=$((OBSLEV_GOOD-OBSLEV))
echo RAISE_H $RAISE_H

echo using disc radius $discR
echo energy $ENERGY
nSamples=100
echo nSamples $nSamples

OUTPUT_FOLDER=/home/enpaudel/icecube/triggerStudy/simFiles/
# DATASETUNIQUE=dataSetUniqueWFRT
DATASETUNIQUE=dataSetUnique
# DATASETUNIQUE=dataSetUnique1_6
# DATASETGEN=dataSetGen1_6
DATASETGEN=dataSetGen
# DATASETUNIQUE=dataSetUniqueFRT
start_time=$SECONDS

#calculate the seed
echo calculating the seed
echo $corsikaRunID
echo $PRIMARY_NAME
echo $pSEED
echo $ENERGY
###########old seed with pSEED 14,5626,402,1608######
# eSEED=$(python -c "import math;print(int(math.log10($ENERGY)))")
# SEED=$(($corsikaRunID+($pSEED*10+$eSEED)*10800000))
#####################################################
eSEED=$(python -c "import math;print(int(math.log10($ENERGY)*10))")
SEED=$(($corsikaRunID+($pSEED*100+$eSEED)*1000000)) #pseed:1,eseed:2,runID:6 digits (total:9 digits)
echo SEED $SEED


DETECTOR_PY=$I3SRC/simprod-scripts/resources/scripts/detector.py
FLAGS2="--UseGSLRNG --gcdfile ${GCD} --noInIce --LowMem --seed ${SEED} --nproc 1 --DetectorName IC86.2019 --no-FilterTrigger"
FLAGS2+=" --inputfile $OUTPUT_FOLDER/${DATASETGEN}/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2"
FLAGS2+=" --output $OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2"

echo running $DETECTOR_PY
echo flags $FLAGS2

$ICETRAY_ENV ${DETECTOR_PY} $FLAGS2
# if [ ! -f $OUTPUT_FOLDER/dataSet/${CORSIKA_ID}GenDet.i3.bz2 ]; then
# 	$ICETRAY_ENV ${DETECTOR_PY} $FLAGS2
# fi
time_det=$(($SECONDS - $start_time))


PHOTONDIR="/cvmfs/icecube.opensciencegrid.org/data/photon-tables"
FILTERING_PY=$I3SRC/filterscripts/resources/scripts/SimulationFiltering.py
FLAGS3=" --needs_wavedeform_spe_corr --photonicsdir ${PHOTONDIR} -g ${GCD}"
FLAGS3+=" --input $OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2"
FLAGS3+=" --output $OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDetFilt.i3.bz2"

echo running $FILTERING_PY
echo flags $FLAGS3

$ICETRAY_ENV ${FILTERING_PY} $FLAGS3
# if [ ! -f $OUTPUT_FOLDER/dataSet/${CORSIKA_ID}GenDetFilt.i3.bz2 ]; then
# 	$ICETRAY_ENV ${FILTERING_PY} $FLAGS3
# fi

time_filt=$(($SECONDS - $start_time))

PROCESS_PY=$I3SRC/filterscripts/resources/scripts/offlineL2/process.py
FLAGS4=" -s --photonicsdir ${PHOTONDIR} -g ${GCD}"
FLAGS4+=" --input $OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDetFilt.i3.bz2"
FLAGS4+=" --output $OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDetFiltProc.i3.bz2"

echo running $PROCESS_PY
echo flags $FLAGS4

echo "$ICETRAY_ENV ${PROCESS_PY} $FLAGS4"
$ICETRAY_ENV ${PROCESS_PY} $FLAGS4

# echo "time_diff=[$SECONDS - $start_time]"
time_proc=$(($SECONDS - $start_time))
# echo $time_diff
# printf "$ENERGY $time_zip $time_gen $time_det $time_filt $time_proc \n" >> $BASEDIR/../energyTimeMulti.txt  
# printf "$ENERGY $time_proc \n" >> $BASEDIR/../energyTimeMultiStudyProper.txt

UNIQUE_PY=$BASEDIR/addUniqueRunID.py
INPUT_FILE=$OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDetFiltProc.i3.bz2

$ICETRAY_ENV $UNIQUE_PY $INPUT_FILE
# rm $ifile
echo removing intermediate files
################################# rm $OUTPUT_FOLDER/dataSet/${PRIMARY_NAME}${CORSIKA_ID}Gen.i3.bz2
rm $OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDet.i3.bz2 $OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDetFilt.i3.bz2
rm $OUTPUT_FOLDER/${DATASETUNIQUE}/${PRIMARY_NAME}${CORSIKA_ID}GenDetFiltProc.i3.bz2

