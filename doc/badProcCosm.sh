#!/bin/bash

: ${PREFDIR:="."}

: ${MYWDIR:="/data/devAlg"}

export Lanes=${Lanes:-3}
export Nthread=${Nthread:-3}
export runOPTL=${runOOPTL:-"-b --run"}
export DetList=${DetList:-"ITS,ITS-TRD,ITS-TOF,ITS-TRD-TOF"}
export TrackSrc=${TrackSrc:-"cosmic-ext"}

: ${Config:=/home/shahoian/alice/O2Align/doc/ITS_RB_all.json}

: ${MAXTF:=-1}

: ${ConfString:="AlignParams.dofConfigJson=$Config;AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=4;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=4;AlignParams.MVTimeSlotsJson=/data/devAlg/MV_B0BPBM.json;$EXTRACONF;"}


script=${MYWDIR}/proc0.sh

echo "Basedir=${MYWDIR}, script: $script"

mkdir -p ${PREFDIR}/milleCosmB0
pushd ${PREFDIR}/milleCosmB0
DataDir="/data/devAlg/data/cosmics_b0/b0_569489r/job_000/" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script #| tee rr1.log
popd

exit
