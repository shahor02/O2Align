#!/bin/bash

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
  env "$@" /data0/devAlg/tests/sol.sh
  popd
}

#export EXTRACONF="AlignParams.statInFile=/data0/devAlg/tests/MilleV8/procStat.root"
#runProc /data0/devAlg/tests/solA2v5bow Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v5bow.json MPSTEER=/data0/devAlg/tests/mpSteer.txt
#runProc /data0/devAlg/tests/solA2v6 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v6.json MPSTEER=/data0/devAlg/tests/mpSteer.txt
#runProc /data0/devAlg/tests/solA2v7 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v7.json MPSTEER=/data0/devAlg/tests/mpSteer.txt

#runProc /data0/devAlg/tests/solA2v8 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v8.json MPSTEER=/data0/devAlg/tests/mpSteerAll.txt
#runProc /data0/devAlg/tests/solA2v8ns Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v8ns.json MPSTEER=/data0/devAlg/tests/mpSteerAll.txt
#runProc /data0/devAlg/tests/solA2v9 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v9.json MPSTEER=/data0/devAlg/tests/mpSteerAll.txt
#runProc /data0/devAlg/tests/solA2v10 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v10.json MPSTEER=/data0/devAlg/tests/mpSteer.txt
#runProc /data0/devAlg/tests/solA2v11 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v11.json MPSTEER=/data0/devAlg/tests/mpSteerV11.txt


#epn000    runProc /data0/devAlg/tests/solA2v12 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v12.json MPSTEER=/data0/devAlg/tests/mpSteerV12.txt

#runProc /data0/devAlg/tests/solA2v13b Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v13b.json MPSTEER=/data0/devAlg/tests/mpSteer.txt
#runProc /data0/devAlg/tests/solA2v13bNO Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v13b.json MPSTEER=/data0/devAlg/tests/mpSteerNoOut.txt
#runProc /data0/devAlg/tests/solA2v13bNOCHI Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v13b.json MPSTEER=/data0/devAlg/tests/mpSteerNoOutChi.txt
#runProc /data0/devAlg/tests/solA2IBv0 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_IBv0.json MPSTEER=/data0/devAlg/tests/mpSteer.txt
#runProc /data0/devAlg/tests/solA2IBv1 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_IBv1.json MPSTEER=/data0/devAlg/tests/mpSteer.txt
#runProc /data0/devAlg/tests/solA2IBv2 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_IBv2.json MPSTEER=/data0/devAlg/tests/mpSteer.txt
#runProc /data0/devAlg/tests/solA2IBv3 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_IBv3.json MPSTEER=/data0/devAlg/tests/mpSteer.txt

#runProc /data0/devAlg/tests/solA2v19s Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v19s.json MPSTEER=/data0/devAlg/tests/mpSteerAll.txt
runProc /data0/devAlg/tests/solA2v19rmc Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v19rmc.json MPSTEER=/data0/devAlg/tests/mpSteerV12.txt

