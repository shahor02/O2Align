#!/bin/bash

export MYWDIR="/data0/devAlg/tests"

export DataList="${MYWDIR}/collCheckList.txt"
export PREFDIR="./"

ErrDet=";GPU_rec_trd.minTrackPt=0.3;;GPU_rec_trd.trkltResRPhiIdeal=0.2;GPU_rec_trd.trkltResVsTanPhiMisalign=1.;AlignParams.extraClsErrYTRD=0.05;AlignParams.extraClsErrZTRD=1.5;AlignParams.extraClsErrYTOF=0.5;AlignParams.extraClsErrZTOF=0.5;"
for ii in $(seq 0 6) ; do
 ErrDet+="AlignParams.extraClsErrYITS[$ii]=0.0003;AlignParams.extraClsErrZITS[$ii]=0.0003;"
done

: ${TRDSHIFT:="AlignParams.applyTRDshift=false"}
export EXTRACONF=${EXTRACONF:-"$TRDSHIFT;AlignParams.useStableRef=false;AlignParams.debugTracksFraction=1.;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=5;AlignParams.minPt=0.8;AlignParams.meanPtB0=0.45;AlignParams.meanPtB0Cosm=4;AlignParams.corrType=2;AlignParams.maxChi2Ndf=15.f;$ErrDet;"}
if [[ -d $(pwd)/GLO/Config/GeometryAligned ]] ; then
    export REMAP="--condition-remap file://$(pwd)=GLO/Config/GeometryAligned"
fi
export OutType="ProcessData"
export MAXTF=${MAXTF:-6000}

${MYWDIR}/badProcTrdTof.sh

