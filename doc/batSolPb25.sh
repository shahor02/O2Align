#!/bin/bash

: ${MYWDIR:="/data0/devAlg/tests"}

runProc() # args: destDir var=value [var=value ...]
{
  if [[ $# -lt 2 ]]; then
    echo "runProc expects arguments: <destDir> <var=value> [<var=value> ...]"
    exit 1
  fi
  local destDir=$1
  shift
  for assignment in "$@"; do
    if [[ ! $assignment =~ ^[a-zA-Z_][a-zA-Z0-9_]*=.*$ ]]; then
      echo "runProc expects variable assignments in the form <var>=<value>: $assignment"
      exit 1
    fi
  done
  echo "Running alignment with assignments $* in dir. $destDir"
  [[ ! -d $destDir ]] && mkdir $destDir
  pushd $destDir
  env "$@" /data0/devAlg/tests/solPb.sh >& procSol.log
  popd
}

export STARTGEOM=${STARTGEOM:-"/data0/devAlg/tests/best/IBEdge/"}
export DetList=${DetList:-"ITS"}
#export SKIPCHECKRESID=1
export MVSLOTS=${MVSLOTS:-"/data1/PbPb/newAlg/MV_B0BPBM.json"}

#DST=/data0/devAlg/tests/solPbHBA0
#runProc $DST PEDEDONE=1 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_HBonly.json MPSTEER=/data0/devAlg/tests/mpSteerPbMP2.txt

runProc /data0/devAlg/tests/solPbHBA0v2  Config=/home/shahoian/alice/O2Align/doc/ITS_RB_HBonly_v2.json MPSTEER=/data0/devAlg/tests/mpSteerPbMP2.txt
runProc /data0/devAlg/tests/solPbHBA0v2b Config=/home/shahoian/alice/O2Align/doc/ITS_RB_HBonly_v2b.json MPSTEER=/data0/devAlg/tests/mpSteerPbMP2.txt

#pushd $DST
#REMAP= /data1/PbPb/batResRS.sh
#/data1/PbPb/reResidAll.sh
#popd
