#!/bin/bash

HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

EXE=$BASEDIR+runDetectorCluster.sh
fileList="$@"
for ifile in $fileList;do
	echo $EXE $ifile
	$EXE $ifile
done