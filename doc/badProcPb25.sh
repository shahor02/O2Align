#!/bin/bash

truncate_log() {
  $O2Align_ROOT/share/macro/truncateLog.sh "$@"
}

echo "/data1/PbPb/newAlg" > "/data1/PbPb/newAlg/list_coll.txt"
: ${DataList:="/data1/PbPb/newAlg/list_coll.txt"}

: ${PREFDIR:="Mille_Pb_MP2"}
: ${MYWDIR:="/data0/devAlg/tests"}
: ${OUTDIR:="IBEdgeS4R"}

export Lanes=15
export Nthreads=3
export runOPTL="-b --run"
#export DetList=${DetList:-"ITS-TRD,ITS-TOF,ITS-TRD-TOF"}
export DetList=${DetList:-"ITS"}
export OutType=${OutType:-"MilleData,MilleSteer"}

: ${Config:=/home/shahoian/alice/O2Align/doc/ITS_RB_all.json}

: ${MAXTF:=-1}

export O2_DPL_MVBIAS=${O2_DPL_MVBIAS:-";mvbias.xyz[0]=0.033;mvbias.xyz[1]=-0.1825;mvbias.xyz[2]=0.0;mvbias.slopeX=-4.27e-05;mvbias.slopeY=-6.96e-05;"}

#ConfString="AlignParams.dofConfigJson=$Config;AlignParams.ITSOverlapMargin=-1;AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=6;AlignParams.MVTimeSlotsJson=/data0/devAlg/data/MV_B0BPBM.json;"
ErrDet=";GPU_rec_trd.minTrackPt=0.3;;GPU_rec_trd.trkltResRPhiIdeal=0.2;GPU_rec_trd.trkltResVsTanPhiMisalign=1.;AlignParams.extraClsErrYTRD=0.05;AlignParams.extraClsErrZTRD=1.5;AlignParams.extraClsErrYTOF=0.5;AlignParams.extraClsErrZTOF=0.5;"
for ii in $(seq 0 6) ; do
 ErrDet+="AlignParams.extraClsErrYITS[$ii]=0.0003;AlignParams.extraClsErrZITS[$ii]=0.0003;"
done

export EXTRACONF=${EXTRACONF:-"AlignParams.MVTimeSlotsJson=/data1/PbPb/newAlg/MV_B0BPBM.json;AlignParams.useStableRef=false;AlignParams.debugTracksFraction=1;AlignParams.maxTracksPerMultiTrackPV=10;AlignParams.usePVConstraintMinTracks=20;AlignParams.useMultiTrackPVConstraint=4;AlignParams.minPt=0.8;AlignParams.meanPtB0=0.45;AlignParams.meanPtB0Cosm=4;AlignParams.corrType=2;AlignParams.maxChi2Ndf=15.f;$ErrDet;"}
export EXTRAOPT=${EXTRAOPT:-}

export REMAP=${REMAP:-"--condition-remap file://${MYWDIR}/best/IBEdge/=GLO/Config/GeometryAligned"}

script=/data0/devAlg/tests/proc0.sh

echo "Basedir=${MYWDIR}, script: $script"

for inps in $(cat $DataList) ; do
  [[ $inps == "#*" ]] && continue
  outs=`echo $inps | sed "s;alien://;;g" | sed "s;\/;_;g"`
  outd="${PREFDIR}/${OUTDIR}/$outs"
  logFile="proc.log"
  maxLogSize="2M"
  echo "Processing $inps to $outd"
  mkdir -p $outd
  pushd $outd
    DataDir="$inps" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script >& $logFile
#    [[ -f mp2data_0.bin ]] && merge mp2data.bin mp2data_*.bin
    [[ -f debug_alg_0.root ]] && merge debug_alg.root debug_alg_*.root
    [[ -f procStat_0.root ]] && merge procStat.root procStat_*.root
    [[ -f procStat.root ]] && truncate_log $logFile 2MB   # truncate only if the job has succeeded.
  popd
done
