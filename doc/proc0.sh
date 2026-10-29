#!/bin/bash

: ${MYWDIR:="/data/devAlg"}

source ~/alice/O2Align/install/init.sh

#Lanes=3 Nthreads=5 runOPTL="-b --run" ./proc0.sh  | tee rr1.log

if [[ -z $DataDir ]] || [[ -z $Config ]] ; then
  echo "Config $Config or DataDir $DataDir are nor defined"
  exit 1  
fi

export ALICEO2_CCDB_CURL_TIMEOUT_DOWNLOAD=200

: ${ConfString:="AlignParams.dofConfigJson=$Config;AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=4;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=6;AlignParams.MVTimeSlotsJson=/scratch/rs/algDev/data/MV_B0BPBM.json;$EXTRACONF;"}
: ${CCDB:=}
: ${Lanes:=1}
: ${TFTHROT:=$Lanes}
: ${DetList:=ITS}
: ${TrackSrc:=ITS,ITS-TPC,ITS-TPC-TRD,ITS-TPC-TRD-TOF,ITS-TPC-TOF}
: ${OutType:="MilleData,MilleSteer"}
: ${EXTRAOPT:=}

echo "CCDB remap: $CCDB"
echo "ConfString: $ConfString"
echo "Lanes: $Lanes"

export DPL_DISABLE_TPC_TRIGGER_READER=1
export DPL_REPORT_PROCESSING=1

if [[ $OutType == *"MilleRes"* ]] || [[ $OutType == "MilleSteer" ]] || [[ $OutType == *"WeakModesDiag"* ]] ; then
    : ${TimeStamp:=1746914681240}
      o2-dev-alignment-workflow ${CCDB} ${EXTRAOPT} --shm-segment-size ${SHMSIZE:-8000000000} --disable-mc  --detectors ${DetList}  --output $OutType --hbfutils-config hbfutils --configKeyValues "HBFUtils.startTime=$TimeStamp;${ConfString};" ${runOPTL}
      echo "Extracted $OutType"
      exit 0
fi

GLOSET="--shm-segment-size ${SHMSIZE:-6000000000} --hbfutils-config ${DataDir}/o2_tfidinfo.root,upstream --timeframes-rate-limit ${TFTHROT} --timeframes-rate-limit-ipcid $(($(date +%s%3N) % 1234567))"

cmd="o2-reader-driver-workflow $GLOSET --max-tf ${MAXTF:--1} | \
o2-dev-alignment-workflow $GLOSET ${CCDB} ${EXTRAOPT} --disable-mc --pipeline barrel-alignment:${Lanes}  --input-dir $DataDir --nthreads ${Nthreads:-1} --output MilleData,MilleSteer \
--track-sources ${TrackSrc} --detectors ${DetList}  --config-macro ${CONFMACRO:-none} \
--configKeyValues \"$PVOPT;$ConfString\" ${runOPTL} "

echo "Executing: $cmd"
echo "$cmd" > wf.log

eval $cmd
# | tee mp2.log

