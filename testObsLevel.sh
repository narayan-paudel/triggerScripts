#!/bin/bash


CORSIKA_FILE_PATH=$1
CORSIKA_FOLDER=$(dirname $(realpath -s $1))
CORSIKA_FOLDER=$CORSIKA_FOLDER/
CORSIKA_ID_BZ=$(basename $1)
CORSIKA_ID=${CORSIKA_ID_BZ%.*}
corsikaRunID=$(echo $CORSIKA_ID | sed 's/[^0-9]*//g'|sed 's/^0*//')
# echo using runID $corsikaRunID
logFile=$CORSIKA_FOLDER$CORSIKA_ID.log
# logFile=/data/sim/IceTop/2009/generated/CORSIKA-ice-top/10410/5.2/DAT051453.log
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
	if (( $OBSLEV > 2834 )); then
		echo $PRIMARY $OBSLEV
	fi
fi

OBSLEV_GOOD=2840
RAISE_H=$((OBSLEV_GOOD-OBSLEV))
echo $RAISE_H

# eSEED=$(python -c "import math;print(int(math.log10($ENERGY)))")
eSEED=$(python -c "import math;print(int(math.log10($ENERGY)*10))")
echo eseed $eSEED
echo corsikaRunID $corsikaRunID
echo pSEED $pSEED
SEED=$(($corsikaRunID+($pSEED*100+$eSEED)*1000000))
# SEED=$(($corsikaRunID+($pSEED*10+$eSEED)*10800000))
# SEED=$(($corsikaRunID+($pSEED*10+$eSEED)*10000000))
# eSEED=$(python -c "import math;print(int(math.log10($ENERGY)*10))")
echo eseed $eSEED
# SEED=$(($corsikaRunID+($pSEED*10+$eSEED)*100000000)) #
echo $SEED