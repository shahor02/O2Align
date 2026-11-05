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

: ${DataList:="/data1/cosm/list_cosm2025_mrg.txt"}

: ${PREFDIR:="."}

: ${MYWDIR:="/data0/devAlg/tests"}

export Lanes=${Lanes:-6}
export Nthread=${Nthread:-5}
export runOPTL=${runOPTL:-"-b --run"}
export DetList=${DetList:-"ITS"} #,ITS-TRD,ITS-TOF,ITS-TRD-TOF"}
export TrackSrc=${TrackSrc:-"cosmic-ext"}

: ${Config:=/home/shahoian/alice/O2Align/doc/ITS_RB_all.json}

: ${MAXTF:=-1}

: ${ConfString:="AlignParams.debugTracksFraction=1.;AlignParams.cosmiGoldOnly=false;AlignParams.dofConfigJson=$Config;AlignParams.maxChi2Ndf=15;AlignParams.minITSClsTrack=3;$EXTRACONF;"}


script=${MYWDIR}/proc0.sh

echo "Basedir=${MYWDIR}, script: $script"

for inps in $(cat $DataList) ; do
  outs=`echo $inps | sed "s;alien://;;g" | sed "s;\/;_;g"`
  outd="${PREFDIR}/cosm/$outs"
  echo "Processing $inps to $outd"
  mkdir -p $outd
  pushd $outd
    DataDir="$inps" Config="$Config" ConfString=${ConfString} MAXTF="$MAXTF" $script #| tee rr1.log
    [[ -f mp2data_0.bin ]] && merge mp2data.bin mp2data_*.bin
    [[ -f debug_alg_0.root ]] && merge debug_alg.root debug_alg_*.root
    [[ -f procStat_0.root ]] && merge procStat.root procStat_*.root
  popd
done

exit
