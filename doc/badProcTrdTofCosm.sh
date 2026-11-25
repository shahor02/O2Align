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

truncate_log() {
  $O2Align_ROOT/share/macro/truncateLog.sh "$@"
}

: ${DataList:="/data1/cosm1/list_cosm_mrg.txt"}

: ${PREFDIR:="Mille_IBEdgeS4R"}

: ${MYWDIR:="/data0/devAlg/tests"}
: ${OUTDIR:="cosmHP4LUT"}

export Lanes=${Lanes:-36}
export Nthread=${Nthread:-5}
export runOPTL=${runOPTL:-"-b --run"}
export DetList=${DetList:-"ITS-TRD,ITS-TOF,ITS-TRD-TOF,TRD-TOF"}
export TrackSrc=${TrackSrc:-"cosmic-ext"}

export REMAP=${REMAP:-"--condition-remap file://${MYWDIR}/best/IBEdge/=GLO/Config/GeometryAligned"}

: ${Config:="/home/shahoian/alice/O2Align/doc/ITSnochipTRDTOF_RB.json"}

: ${MAXTF:=-1}

ErrDet=";GPU_rec_trd.minTrackPt=0.3;;GPU_rec_trd.trkltResRPhiIdeal=0.2;GPU_rec_trd.trkltResVsTanPhiMisalign=1.;AlignParams.extraClsErrYTRD=0.05;AlignParams.extraClsErrZTRD=1.5;AlignParams.extraClsErrYTOF=0.5;AlignParams.extraClsErrZTOF=0.5;"
for ii in $(seq 0 6) ; do
 ErrDet+="AlignParams.extraClsErrYITS[$ii]=0.0003;AlignParams.extraClsErrZITS[$ii]=0.0003;"
done

export EXTRACONF=${EXTRACONF:-"AlignParams.debugTracksFraction=0.1;AlignParams.maxTracksPerMultiTrackPV=4;AlignParams.minPt=0.6;AlignParams.meanPtB0=0.45;AlignParams.meanPtB0Cosm=4;AlignParams.corrType=2;AlignParams.maxChi2Ndf=14.f;$ErrDet;"}
export EXTRAOPT=${EXTRAOPT:-}

: ${ConfString:="AlignParams.cosmiGoldOnly=false;AlignParams.dofConfigJson=$Config;AlignParams.minITSClsTrack=3;$EXTRACONF;"}


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
    [[ -f mp2data_0.bin ]] && merge mp2data.bin mp2data_*.bin
    [[ -f debug_alg_0.root ]] && merge debug_alg.root debug_alg_*.root
    [[ -f procStat_0.root ]] && merge procStat.root procStat_*.root
    [[ -f procStat.root ]] && truncate_log $logFile 2MB   # truncate only if the job has succeeded.
  popd
done

exit
