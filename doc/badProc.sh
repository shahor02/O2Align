#!/bin/bash

: ${PREFDIR:="."}

export Lanes=15
export Nthreads=3
export runOPTL="-b --run"

: ${Config:=/home/shahoian/alice/O2Align/doc/ITS_RB_all.json}

: ${MAXTF:=-1}

#ConfString="AlignParams.dofConfigJson=$Config;AlignParams.ITSOverlapMargin=-1;AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=6;AlignParams.MVTimeSlotsJson=/data0/devAlg/data/MV_B0BPBM.json;"


script=/data0/devAlg/tests/proc0.sh

mkdir -p ${PREFDIR}/milleB0
pushd ${PREFDIR}/milleB0
DataDir="/data0/devAlg/data/B0" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script  | tee rr1.log
popd

mkdir -p ${PREFDIR}/milleBPlus
pushd ${PREFDIR}/milleBPlus
DataDir="/data0/devAlg/data/BPlus" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script  | tee rr1.log
popd

mkdir -p ${PREFDIR}/milleBMinus
pushd ${PREFDIR}/milleBMinus
DataDir="/data0/devAlg/data/BMinus" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script  | tee rr1.log
popd
