#!/bin/bash


for j in $(seq 122704204.0 122704306.0);do
	echo condor_rm $j
	condor_rm $j;
done

# condor_rm 