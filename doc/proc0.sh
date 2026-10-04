#!/bin/bash

#Lanes=3 Nthreads=5 runOPTL="-b --run" ./proc0.sh  | tee rr1.log

: ${DataDir:="data/_alice_data_2025_LHC25ab_562862_apass1_alignment_its_0010_o2_ctf_run00562862_orbit0261432416_tf0000451373_epn313"}

: ${confString:="AlignParams.dofConfigJson=/home/shahoian/alice/O2Align/doc/ITS_RB_from_staves.json"}
#export ALICEO2_CCDB_LOCALCACHE=ccdb

export DPL_DISABLE_TPC_TRIGGER_READER=1
export DPL_REPORT_PROCESSING=1
GLOSET="--shm-segment-size ${SHMSIZE:-8000000000} --hbfutils-config ${DataDir}/o2_tfidinfo.root,upstream --timeframes-rate-limit ${TFTHROT:-1} --timeframes-rate-limit-ipcid $(($(date +%s%3N) % 1234567))"

: ${PVOPT:="pvertexer.useMeanVertexConstraint=false;pvertexer.meanVertexExtraErrSelection=0.2;pvertexer.iniScale2=100;pvertexer.acceptableScale2=10."}

cmd="o2-reader-driver-workflow $GLOSET --max-tf ${MAXTF:--1} | \
o2-dev-alignment-workflow $GLOSET --disable-mc --pipeline barrel-alignment:${Lanes:-1}  --input-dir $DataDir --nthreads ${Nthreads:-1} --output MilleData,MilleSteer \
--track-sources ${TrackSrc:-ITS,ITS-TPC,ITS-TPC-TRD,ITS-TPC-TRD-TOF,ITS-TPC-TOF} --detectors ${DetList:-ITS}  --config-macro ${CONFMACRO:-none} \
--configKeyValues \"$PVOPT;$confString\" ${runOPTL} "

echo "Executing: $cmd"
echo "$cmd" > wf.log

eval $cmd
# | tee mp2.log

