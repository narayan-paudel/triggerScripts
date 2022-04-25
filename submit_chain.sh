#!/bin/env bash

function print_help {
    cat <<EOF 1>&2
 usage: ${0} <output> <input>... [ -0 DETOPTIONS -1 LV1OPTIONS -2 LV2OPTIONS -3 LV3OPTIONS -ac -d DETECTOR -f NFAMRES -g GCD -h -i INITIALRUN -j JUMPFILE -m MCDATASET -n NSAMPLES -o ITSOPTIONS -r RUNS -s DATASET -t ]
        
   commands:
     <output>       Output folder name. The output file(s) will be <output>/generates/(topsimulator|detector)/[TopSimulator|Detector]_DETECTOR_corsika_icetop.MCDATASET.RUN.i3.bz2
                    for level 0s, <output>/filtered/Level[1|2]_DETECTOR_corsika_icetop.MCDATASET.RUN.i3.bz2, for level 1 and 2, Level3_DETECTOR_DATASET_Run.RUN.i3.bz2 for level3
                    The condor files will be saved in the realtive condor folders
     <input>        Folder where find CORISKA inputs (DATxxxxxx format); ideally is expected an energy folder. i.e. <inputpath>/6.0
     
   options:
     -0 DETOPTIONS  Enquoted other options for detector.py (level 0 simulation)
     -1 LV1OPTIONS  Enquoted other options for SimulationFiltering.py (level 1 simulation)
     -2 LV2OPTIONS  Enquoted other options for process.py (level 2 simulation)
     -3 LV3OPTIONS  Enquoted other options for level3_iceprod.py (level 3 simulation); available only if -a is declared
     -a             Analysis level, i.e. include level 3 in the production
     -c             Cut data with the trigger filtering. Useful for Level 0, 1 and 2
     -d DETECTOR    Detector type/geometry (IC79, IC86, IC86.2012...): the full name will be used only for L3 reconsteuction,
                     while for the other only the base (i.e. IC86 if declared IC86.2012) [default: IC86.2012]
     -f NFAMRES     How many frames will be processed (0=unbound) [default: 0]
     -g GCD         GCD file [default: /data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz]
     -h             Print this help
     -i INITIALRUN  Set initial run of the CORSIKA to start the processes (see '-r' to end after given counts) [default: 1]
     -j JUMPFILE    File (.txt) containing input file to NOT BE submitted.
     -m MCDATASET   Number of corsika dataset [default: 10410]
     -n NSAMPLES    How many samples per run will be simulated by icetopshowergenerator.py [default: 100]
     -o ITSOPTIONS  Enquoted other options for icetopsimulator.py; example: -o "--raise-observation-level 3"
     -r RUNS        Set maximum run id (this included) to be processed. -1 = all runs [default: -1]
     -s DATASET     Number of dataset [default: 12360]
     -t             Test: launch everything as a DryRun

${1}
EOF
    exit ${2:0}
}

print_help

# function make_condor_steering {
#     condor=${1-test}
#     out=${2-test}
#     args=${3-}
#     cat <<EOF > ${condor}.sub

# #############################################
# #this shell script will call my python script for condor 
# #  (condor is the batch scheduler for NPX3)
# # This takes variables defined in the DAG file whose names 
# #   are case-insentivite. (here: infile, outfile, jobname, gcd, datatype)
# # Run using command $condor_submit condor.submit
# # (if necess, kill with $ condor_rm jobnum
# #      or $ condor_rm myusername
# ###############################################
# # Jobname will name your output data, .out, .log, .err files
# # This script name will show up on the cluster queue:

# executable = ${condor}.sh

# ####################################
# # Once you've set the stuff below here to your directories, 
# # it probably won't change from job to job
# ###################################

# #this is where the log, out, err files will live

# output         = ${condor}.out
# error          = ${condor}.err
# log            = ${out}.log

# # This option takes all your env variables from your current session
# # so there's no need to submit your env-shell.sh to the cluster

# getenv         = true
# universe       = vanilla
# notification   = never

# # Condor 7.8 (new on npx4 Dec 2012) comes with a better implementation of "dynamic slots".  
# #request_cpus = <num_cpus>

# request_memory = 32GB

# #request_disk = <disk_in_KB>
# # to the condor submit file.  A dynamic slot is generated with the requested resources.  If your condor submit file is created without these new request lines a generic slot with 1GB of RAM and 1CPU core is generated similar to todays setup.

# ########################(drum roll..... SUBMIT!)####################

# Arguments = ${args}
# queue 1

# EOF
# }

# function make_condor_exec {
#     condor=${1-test}
#     cmdl=${2-"echo ./test"}
# #    storename=${3-/scratch/${USER}/test}
# #    storedir=$(dirname ${storename})
#     destpath=$(dirname ${condor})
#     exename=(${cmdline})

#     cat<<EOF > ${condor}.sh
# #!/bin/env bash

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PWD}

# cd ${PWD}

# echo "- '${exename[1]}' starts"
# ${cmdl}
# echo " => ${exename[1]} returned" \$?

# EOF
#     chmod +x ${condor}.sh
# }

# function make_condor_chain {
#     condor=${1-test}
#     cmdlITS=${2}
#     cmdlDet=${3}
#     cmdlLv1=${4}
#     cmdlLv2=${5}
#     cmdlLv3=${6-""}
#     datinput=${7-"./DAT000001"}
#     if [ "${cmdlLv3}" == "" ]
#     then
# 	doreco=0
#     else
# 	doreco=1
#     fi
#     destpath=$(dirname ${condor})
#     datbase=`basename ${datinput}`
#     datext=${datbase##*.}
#     datfile=${datbase%.*}
#     exename=(${cmdline})
#     cat<<EOF > ${condor}.sh
# #!/bin/env bash

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PWD}

# cd ${PWD}

# resC=1
# resD=1
# res1=1
# res2=1

# echo "- Check and copy CORSIKA file"
# if [ ! -e "${destpath}/${datbase}" ] && [ ! -e "${destpath}/${datfile}" ]
# then
#     cp ${datinput} ${destpath}
# fi
# if [ "${datext}" == "bz2" ]
# then
#     if [ ! -e "${destpath}/${datfile}" ]
#     then
# 	bzip2 -d ${destpath}/${datfile}.${datext}
#     fi
#     if [ -e "${destpath}/${datfile}.${datext}" ]
#     then
#         rm ${destpath}/${datfile}.${datext}
#     fi
# fi

# echo "- Chain starts"
# ${cmdlITS}
# resC=\$?
# if [ \${resC} == 0 ]
# then
#     rm ${destpath}/${datfile}
#     ${cmdlDet}
#     resD=\$?
# else
#     echo " => TopSimulator returned" \$? 
# fi
# if [ \${resD} == 0 ]
# then
#     ${cmdlLv1}
#     res1=\$?
# else
#     echo " => Detector returned" \$? 
# fi
# if [ \${res1} == 0 ]
# then
#     ${cmdlLv2}
#     res2=\$?
# else
#     echo " => Level1 returned" \$? 
# fi
# if [ ${doreco} == 1 ]
# then
#     if [ \${res2} == 0 ]
#     then
#         ${cmdlLv3}
#         echo " => Chain (up to Level 3) returned" \$?
#     else
#         echo " => Level2 returned" \${res2}
#     fi
# else
#     echo " => Chain (up to Level 2) returned" \${res2}
# fi

# EOF
#     chmod +x ${condor}.sh
# }

# function get_ebin {
#     path=${1-5.0}
#     ebin=`basename ${path}`
#     ebin=${ebin%.*}
#     echo $((ebin-5))
# }

# function get_radius {
#     ebin=${1-0}
#     base=${2-800}
#     echo $((base + ebin*300 + (2*ebin/3)*300 + (ebin/3)*300))
# }

# function get_i3_runid {
#     path=${1-test.000001.i3.bz2}
#     name=`basename ${path}`
#     name=${name%.i3*}
#     num=${name##*.}
#     #name=${name//[!0-9]/}
#     echo ${num} | sed 's/^0*//'
# }

# function get_corsika_runid {
#     path=${1-DAT000001}
#     name=`basename ${path}`
#     name=${name%%.*}
#     name=${name//[!0-9]/}
#     echo ${name} | sed 's/^0*//'
# }


# scratch=/scratch/${USER}/icetopsim
# basepath=${PWD}

# # default args
# DATASET=12360
# GCDdefault="/data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_2012.56063_V1_OctSnow.i3.gz"
# GCD=${GCDdefault}
# NSAMPLES=100
# NFARMES=0
# ITSOPTIONS=""
# DETOPTIONS=""
# LV1OPTIONS=""
# LV2OPTIONS=""
# LV3OPTIONS=""
# INITIALRUN=1
# RUNS=0
# DRYRUN=""
# MCDATASET=10410
# DETECTOR="IC86.2012"
# CUT=0
# PHOTONDIR="/cvmfs/icecube.opensciencegrid.org/data/photon-tables"
# RECO=0
# JUMPFILE=""

# # get args
# while getopts ":0:1:2:3:acd:f:g:hi:j:m:n:o:r:s:t" flag
# do
#     case "$flag" in
# 	0) DETOPTIONS=$OPTARG;;
# 	1) LV1OPTIONS=$OPTARG;;
# 	2) LV2OPTIONS=$OPTARG;;
# 	3) LV3OPTIONS=$OPTARG;;
# 	a) RECO=1;;
# 	c) CUT=1;;
# 	d) DETECTOR=$OPTARG;;
# 	f) NFARMES=$OPTARG;;
# 	g) GCD=$OPTARG;;
#         h) print_help;;
# 	i) INITIALRUN=$OPTARG;;
# 	j) JUMPFILE=$OPTARG;;
# 	m) MCDATASET=$OPTARG;;
# 	n) NSAMPLE=$OPTARG;;
# 	o) ITSOPTIONS=$OPTARG;;
#         r) RUNS=$OPTARG;;
# 	s) DATASET=$OPTARG;;
# 	t) DRYRUN="dryrun!";;
#     esac
# done

# # get positionals
# if (( $#<$OPTIND+1 ))
# then
#     print_help "At least 2 arguments (beyond the options) are required!" 1
# fi
# OUTPUT=${PWD}/${@:$OPTIND:1}
# INPUTS=${@:$OPTIND+1}

# L3DETECTOR=${DETECTOR}
# DETECTOR=${DETECTOR%.*}

# # = adjust parameters based on the LEVEL set
# # calculate nproc
# nproc=0
# for inp in ${INPUTS}
# do
#     for dat in `ls ${inp}/DAT*[0-9]`
#     do
# 	runid=`get_corsika_runid ${dat}`
# 	if [ ${INITIALRUN} -gt ${runid} ]
# 	then
# 	    continue
# 	fi
# 	let nproc+=1
#     done
#     #let nproc+=`ls ${inp}/DAT*[0-9] | wc -l`
# done

# if [ ${RUNS} -gt 0 ] && [ ${nproc} -gt ${RUNS} ]
# then
#     nproc=${RUNS}
# fi

# # calculate seed
# seed=$((1234+DATASET*3))

# # = build structures
# # output
# outname=$(basename $OUTPUT)
# destcondor=${PWD}/${outname}/condor
# mkdir -p ${destcondor}
# # log folder
# tmpdir="${scratch}/${outname}_chain"
# mkdir -p ${tmpdir}

# dset=${DATASET}
# mcset=${MCDATASET}
# DATASET=`printf "%06d" ${DATASET}`
# MCDATASET=`printf "%06d" ${MCDATASET}`

# outputITS="${OUTPUT}/generated/topsimulator"
# mkdir -p ${outputITS}/condor
# outputITS+="/TopSimulator_${DETECTOR}_corsika_icetop.${MCDATASET}."

# outputDet="${OUTPUT}/generated/detector"
# mkdir -p ${outputDet}/condor
# outputDet+="/Detector_${DETECTOR}_corsika_icetop.${MCDATASET}."

# outputLv1="${OUTPUT}/filtered/level1"
# mkdir -p ${outputLv1}/condor
# outputLv1+="/Level1_${DETECTOR}_corsika_icetop.${MCDATASET}."

# outputLv2="${OUTPUT}/filtered/level2"
# mkdir -p ${outputLv2}/condor
# outputLv2+="/Level2_${DETECTOR}_corsika_icetop.${MCDATASET}."

# if [ ${RECO} == 1 ]
# then
#     outputLv3="${OUTPUT}/filtered/level3"
#     mkdir -p ${outputLv3}/condor
#     outputLv3+="/Level3_${L3DETECTOR}_${DATASET}_Run"
# fi

# # set execution in case of a dryrun
# if [ "${DRYRUN}" == "" ]
# then
#     cmdbase="time"
# else
#     cmdbase="echo"
# fi

# # command lines
# cmdlineITS="${cmdbase} ${I3_BUILD}/simprod-scripts/resources/scripts/icetopshowergenerator.py --UseGSLRNG --gcdfile ${GCD} --seed ${seed} --nproc ${nproc} --samples ${NSAMPLES}"
# if [ "${ITSOPTIONS}" != "" ]
# then
#     cmdlineITS+=" ${ITSOPTIONS}"
# fi

# cmdlineDet="${cmdbase} ${I3_BUILD}/simprod-scripts/resources/scripts/detector.py --UseGSLRNG --gcdfile ${GCD} --noInIce --LowMem --seed ${seed} --nproc ${nproc} --DetectorName ${DETECTOR}"
# if [ ${CUT} -eq 0 ]
# then
#     cmdlineDet+=" --no-FilterTrigger"
# fi
# if [ "${DETOPTIONS}" != "" ]
# then
#     cmdlineDet+=" ${DETOPTIONS}"
# fi

# cmdlineLv1="${cmdbase} ${I3_BUILD}/filterscripts/resources/scripts/SimulationFiltering.py --needs_wavedeform_spe_corr --photonicsdir ${PHOTONDIR} -g ${GCD}"
# if [ "${LV1OPTIONS}" != "" ]
# then
#     cmdlineLv1+=" ${LV1OPTIONS}"
# fi

# cmdlineLv2="${cmdbase} ${I3_BUILD}/filterscripts/resources/scripts/offlineL2/process.py -s --photonicsdir ${PHOTONDIR} -g ${GCD}"
# if [ "${LV2OPTIONS}" != "" ]
# then
#     cmdlineLv2+=" ${LV2OPTIONS}"
# fi

# if [ ${RECO} == 1 ]
# then
#     cmdlineLv3="${cmdbase} ${I3_BUILD}/icetop_Level3_scripts/resources/scripts/level3_iceprod.py -m --waveforms --spe-corr --dataset=${dset} -d ${L3DETECTOR}"
#     if [ "${GCD}" != "${GCDdefault}" ]
#     then
# 	gcdbase=`basename ${GCD}`
# 	gcdname=${gcdbase%.i3*}
# 	gcdext=${gcdbase#${gcdname}}
# 	l3gcd=${destcondor}/${gcdname}_L3${gcdext}
# 	if [ ! -e "${l3gcd}" ]
# 	then
# 	    python ${I3_BUILD}/icetop_Level3_scripts/resources/scripts/MakeL3GCD_MC.py --MCgcd ${GCD} --output ${l3gcd}
# 	fi
# 	cmdlineLv3+=" --L2-gcdfile=${GCD} --L3-gcdfile=${l3gcd}"
#     fi
#     if [ "${LV3OPTIONS}" != "" ]
#     then
# 	cmdlineLv3+=" ${LV3OPTIONS}"
#     fi
# else
#     cmdlineLv3=""
# fi

# # events cap
# if [ ${NFARMES} -gt 0 ]
# then
#     cmdlineDet+=" -n ${NFARMES}"    
#     cmdlineLv1+=" -n ${NFARMES}"    
#     cmdlineLv2+=" -n ${NFARMES}"
#     if [ ${RECO} == 1 ]
#     then
# 	cmdlineLv3+=" -n ${NFARMES}"
#     fi
# fi

# if [ "${JUMPFILE}" != "" ] && [ -e ${JUMPFILE} ]
# then
#     jumpers=(`cat  ${JUMPFILE}`)
#     nj=${#jumpers[@]}
#     j=0
#     let nj=nj-1
# fi


# # = RUN!
# echo "= $outname:"
# rm -f ${destcondor}/submit.log

# procnum=1
# for inp in ${INPUTS}
# do
#     if [ "${inp}" == "condor" ]
#     then
# 	continue
#     fi

#     for dat in `ls ${inp}/DAT*`
#     do
# 	iname=`basename ${dat}`
# 	cname=${iname%.*}
# 	ext=${iname##*.}
# 	runid=`get_corsika_runid ${dat}`
# 	if [ "${runid}" == "" ]
# 	then 
# 	    runid=0
# 	fi
# 	if [ ${RUNS} != 0 ] && [ ${procnum} -gt ${RUNS} ]
# 	then
# 	    break
# 	elif [ ${INITIALRUN} -gt ${runid} ]
# 	then 
# 	    continue
# 	elif [ "${ext}" == "log" ] || [ "${ext}" == "long" ] || [ "${ext}" == "lst" ] || [ "${dat}" == "condor" ]
# 	then
# 	    continue
# 	fi
# 	if [ "${JUMPFILE}" != "" ] && [ -e ${JUMPFILE} ]
# 	then
# 	    jname=`basename ${jumpers[j]}`
# 	    jname=${jname%.*}
# 	    if [ "${jname}" == "${cname}" ]
# 	    then
# 		let j+=1
# 		let procnum+=1
# 		continue
# 	    fi
# 	fi
# 	runname=`printf "%06d" ${runid}`
# 	#let runid=${runid}-1

# 	#thisdest="${destcondor}/Run${runname}_${iname}_chain"
# 	thisdest="${destcondor}/Run${runname}"
# 	thiscondor="${tmpdir}/Run${runname}"
# 	make_condor_steering "${thisdest}" "${thiscondor}"

# 	ITSname="${outputITS}${runname}"
# 	ITSdest="$(dirname ${ITSname})/condor/$(basename ${ITSname})"
# 	Detname="${outputDet}${runname}"
# 	Detdest="$(dirname ${Detname})/condor/$(basename ${Detname})"
# 	Lv1name="${outputLv1}${runname}"
# 	Lv1dest="$(dirname ${Lv1name})/condor/$(basename ${Lv1name})"
# 	Lv2name="${outputLv2}${runname}"
# 	Lv2dest="$(dirname ${Lv2name})/condor/$(basename ${Lv2name})"
# 	ITS="${cmdlineITS} --RunID ${runid} --procnum ${procnum} --inputfilelist ${destcondor}/${cname} --outputfile ${ITSname}.i3.bz2 > ${ITSdest}.out 2> ${ITSdest}.err"
# 	Det="${cmdlineDet} --RunID ${runid} --procnum ${procnum} --inputfile ${ITSname}.i3.bz2 --outputfile ${Detname}.i3.bz2 > ${Detdest}.out 2> ${Detdest}.err"
# 	Lv1="${cmdlineLv1} -i ${Detname}.i3.bz2 -o ${Lv1name}.i3.bz2 > ${Lv1dest}.out 2> ${Lv1dest}.err"
# 	Lv2="${cmdlineLv2} -i ${Lv1name}.i3.bz2 -o ${Lv2name}.i3.bz2 > ${Lv2dest}.out 2> ${Lv2dest}.err"
# 	if [ ${RECO} == 1 ]
# 	then
# 	    Lv3name="${outputLv3}${runname}"
# 	    Lv3dest="$(dirname ${Lv3name})/condor/$(basename ${Lv3name})"
# 	    Lv3="${cmdlineLv3} --run=${runid} -i ${Lv2name}.i3.bz2 -o ${Lv3name}.i3.bz2 > ${Lv3dest}.out 2> ${Lv3dest}.err"
# 	fi
# 	make_condor_chain "${thisdest}" "${ITS}" "${Det}" "${Lv1}" "${Lv2}" "${Lv3}" "${dat}"

# 	condor_submit ${thisdest}.sub >> ${destcondor}/submit.log
# 	if [[ $? == 0 ]]
#         then	
# 	    let procnum+=1
# 	fi
#     done    
# done
# counts=$((procnum-1))

# # = summary printout of the submiting
# if [[ ${counts} -eq 1 ]]
# then
#     echo " -> 1 job successfully submitted"
# elif [[ ${counts} -gt 1 ]]
# then
#     echo " -> ${counts} jobs successfully submitted"
# else
#     echo " -> No jobs are successfully submitted"
# fi

# ##### -o "--raise-observation-level 3" for the 2012 dataset
