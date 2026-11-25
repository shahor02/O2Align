#!/bin/bash

merge() # arge: dest_file inp1 [inp2 ...]
{
  if [[ $# -lt 2 ]]; then
    echo "merge expects at least one output and one input file"
    exit 1
  fi
  local outf=$1  
  shift
  if [[ "$outf" == *.root ]] ; then
    echo "executing: hadd $outf $@"
    hadd $outf "$@"
  elif [[ "$outf" == *.bin ]] ; then
    echo "executing: cat $@ > $outf"
    cat "$@" > $outf
  else
    echo "Unknown file types $@ -> $outf provided for merging"
    exit 1  
  fi
  [[ -f $outf ]] && rm "$@"
}

: ${MYWDIR:="/data0/devAlg/tests"}

source ~/alice/O2Align/install/init.sh
#Lanes=3 Nthreads=5 runOPTL="-b --run" ./proc0.sh  | tee rr1.log

if [[ -z $DataDir ]] || [[ -z $Config ]] ; then
  echo "Config $Config or DataDir $DataDir are nor defined"
  exit 1  
fi

echo "REMAP= $REMAP"

export ALICEO2_CCDB_CURL_TIMEOUT_DOWNLOAD=200

: ${ConfString:="AlignParams.dofConfigJson=$Config;${STARTMILLE:+AlignParams.startFromMilleRes=$STARTMILLE;}AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=6;${MVSLOTS:+AlignParams.MVTimeSlotsJson=$MVSLOTS;}AlignParams.statInFile=$MYWDIR/Mille_IBEdgeS4R/procStatColl.root;$EXTRACONF;"}
#: ${ConfString:="AlignParams.dofConfigJson=$Config;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=6;AlignParams.MVTimeSlotsJson=/data0/devAlg/data/MV_B0BPBM.json"}
#: ${ConfString:="AlignParams.dofConfigJson=$Config;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20;"}
#export ALICEO2_CCDB_LOCALCACHE=ccdb

: ${Lanes:=1}
: ${TFTHROT:=$Lanes}
: ${DetList:=ITS}
: ${TrackSrc:=ITS,ITS-TPC,ITS-TPC-TRD,ITS-TPC-TRD-TOF,ITS-TPC-TOF}
: ${OutType:="MilleData,MilleSteer"}

export DPL_DISABLE_TPC_TRIGGER_READER=1
export DPL_REPORT_PROCESSING=1

if [[ $OutType == *"MilleRes"* ]] || [[ $OutType == "MilleSteer" ]] || [[ $OutType == *"WeakModesDiag"* ]] || [[ $OutType == "Regularize" ]]  ; then
    : ${TimeStamp:=1746914681240}
      cmd="o2-dev-alignment-workflow ${REMAP} ${EXTRAOPT} --shm-segment-size ${SHMSIZE:-8000000000} --disable-mc  --detectors ${DetList}  --output $OutType --hbfutils-config hbfutils --configKeyValues \"HBFUtils.startTime=$TimeStamp;${ConfString};\" ${runOPTL}"
      echo "Extracting $OutType with"
      echo "Executing $cmd"
      eval $cmd
      exit 0
fi

GLOSET="--shm-segment-size ${SHMSIZE:-16000000000} --hbfutils-config ${DataDir}/o2_tfidinfo.root,upstream --timeframes-rate-limit ${TFTHROT} --timeframes-rate-limit-ipcid $(($(date +%s%3N) % 1234567))"

cmd="o2-reader-driver-workflow $GLOSET --max-tf ${MAXTF:--1} | \
o2-dev-alignment-workflow ${REMAP} ${EXTRAOPT} $GLOSET --disable-mc --pipeline barrel-alignment:${Lanes}  --input-dir $DataDir --nthreads ${Nthreads:-1} --output ${OutType} \
--track-sources ${TrackSrc} --detectors ${DetList}  --config-macro ${CONFMACRO:-none} \
--configKeyValues \"$PVOPT;$ConfString\" ${runOPTL} "

echo "Executing: $cmd"
echo "$cmd" > wf.log

eval $cmd
[[ -f debug_alg_0.root ]] && merge debug_alg.root debug_alg_*.root
[[ -f procStat_0.root ]] && merge procStat.root procStat_*.root

# | tee mp2.log

