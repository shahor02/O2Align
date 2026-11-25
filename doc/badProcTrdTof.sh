#!/bin/bash

truncate_log() {
  $O2Align_ROOT/share/macro/truncateLog.sh "$@"
}

: ${DataList:="/data0/devAlg/data/list_coll.txt"}

: ${PREFDIR:="Mille_IBEdgeS4R"}
: ${MYWDIR:="/data0/devAlg/tests"}
: ${OUTDIR:="coll_LutPt08ITSSEED"}

export Lanes=${Lanes:-20}
export Nthreads=${Nthreads:-4}
export runOPTL=${runOPTL:-"-b --run"}
#export DetList=${DetList:-"ITS"}
export DetList=${DetList:-"ITS-TRD,ITS-TOF,ITS-TRD-TOF"}
export OutType=${OutType:-"MilleData,MilleSteer"}

: ${Config:=/home/shahoian/alice/O2Align/doc/ITSnochipTRDTOF_RB.json}

: ${MAXTF:=-1}

export O2_DPL_MVBIAS=${O2_DPL_MVBIAS:-";mvbias.xyz[0]=0.033;mvbias.xyz[1]=-0.1825;mvbias.xyz[2]=0.0;mvbias.slopeX=-4.27e-05;mvbias.slopeY=-6.96e-05;"}

#ConfString="AlignParams.dofConfigJson=$Config;AlignParams.ITSOverlapMargin=-1;AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=6;AlignParams.MVTimeSlotsJson=/scratch/rs/algDev/data/MV_B0BPBM.json;"

ErrDet=";GPU_rec_trd.minTrackPt=0.3;;GPU_rec_trd.trkltResRPhiIdeal=0.2;GPU_rec_trd.trkltResVsTanPhiMisalign=1.;AlignParams.extraClsErrYTRD=0.05;AlignParams.extraClsErrZTRD=1.5;AlignParams.extraClsErrYTOF=0.5;AlignParams.extraClsErrZTOF=0.5;"
for ii in $(seq 0 6) ; do
 ErrDet+="AlignParams.extraClsErrYITS[$ii]=0.0003;AlignParams.extraClsErrZITS[$ii]=0.0003;"
done

export EXTRACONF=${EXTRACONF:-"AlignParams.useStableRef=false;AlignParams.debugTracksFraction=1;AlignParams.maxTracksPerMultiTrackPV=5;AlignParams.useMultiTrackPVConstraint=2;AlignParams.minPt=0.8;AlignParams.meanPtB0=0.45;AlignParams.meanPtB0Cosm=4;AlignParams.corrType=2;AlignParams.maxChi2Ndf=15.f;$ErrDet;"}
export EXTRAOPT=${EXTRAOPT:-}
#export EXTRAOPT=${EXTRAOPT:-"--stacktrace-on-signal none"}

export REMAP=${REMAP:-"--condition-remap file://${MYWDIR}/best/IBEdge/=GLO/Config/GeometryAligned"}

script=${MYWDIR}/proc0.sh

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
