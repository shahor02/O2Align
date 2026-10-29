#!/bin/bash

export O2_DPL_MVBIAS=";mvbias.xyz[0]=0.0332;mvbias.xyz[1]=-0.1826;mvbias.xyz[2]=0.0;mvbias.slopeX=-4.27e-05;mvbias.slopeY=-6.96e-05;"

: ${PREFDIR:="."}

: ${MYWDIR:="/data/devAlg"}

export Lanes=${Lanes:-3}
export Nthread=${Nthread:-3}
export runOPTL=${runOOPTL:-"-b --run"}
export DetList=${DetList:-"ITS,ITS-TRD,ITS-TOF,ITS-TRD-TOF"}

: ${Config:=/home/shahoian/alice/O2Align/doc/ITS_RB_all.json}

: ${MAXTF:=-1}

: ${ConfString:="AlignParams.dofConfigJson=$Config;AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=4;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=4;AlignParams.MVTimeSlotsJson=/data/devAlg/MV_B0BPBM.json;$EXTRACONF;"}


script=${MYWDIR}/proc0.sh

echo "Basedir=${MYWDIR}, script: $script"

mkdir -p ${PREFDIR}/milleB0
pushd ${PREFDIR}/milleB0
DataDir="/data/devAlg/data/_alice_data_2025_LHC25ab_562862_apass1_alignment_its_0010_o2_ctf_run00562862_orbit0261432416_tf0000451373_epn313/" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script #| tee rr1.log
popd

exit

mkdir -p ${PREFDIR}/milleB0bis
pushd ${PREFDIR}/milleB0bis
DataDir="/scratch/rs/algDev/data/B0bis" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script  | tee rr1.log
popd

exit 0

mkdir -p ${PREFDIR}/milleBPlus
pushd ${PREFDIR}/milleBPlus
DataDir="/scratch/rs/algDev/data/BPlus" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script  | tee rr1.log
popd

mkdir -p ${PREFDIR}/milleBMinus
pushd ${PREFDIR}/milleBMinus
DataDir="/scratch/rs/algDev/data/BMinus" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script  | tee rr1.log
popd
