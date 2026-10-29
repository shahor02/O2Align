#!/bin/bash

: ${MYWDIR:=/scratch/rs/algDev/tests}

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
  env "$@" ${MYWDIR}/sol.sh
  popd
}

#export EXTRACONF="AlignParams.statInFile=${MYWDIR}/MilleV8/procStat.root"
#runProc ${MYWDIR}/solA2v5bow Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v5bow.json MPSTEER=${MYWDIR}/mpSteer.txt
#runProc ${MYWDIR}/solA2v6 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v6.json MPSTEER=${MYWDIR}/mpSteer.txt
#runProc ${MYWDIR}/solA2v7 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v7.json MPSTEER=${MYWDIR}/mpSteer.txt

#runProc ${MYWDIR}/solA2v8 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v8.json MPSTEER=${MYWDIR}/mpSteerAll.txt
#runProc ${MYWDIR}/solA2v8ns Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v8ns.json MPSTEER=${MYWDIR}/mpSteerAll.txt
#runProc ${MYWDIR}/solA2v9 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v9.json MPSTEER=${MYWDIR}/mpSteerAll.txt
#runProc ${MYWDIR}/solA2v10 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v10.json MPSTEER=${MYWDIR}/mpSteer.txt
#runProc ${MYWDIR}/solA2v11 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v11.json MPSTEER=${MYWDIR}/mpSteerV11.txt


#runProc ${MYWDIR}/solA2v11re Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v11.json MPSTEER=${MYWDIR}/mpSteerV12.txt
#runProc ${MYWDIR}/solA2v12 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v12.json MPSTEER=${MYWDIR}/mpSteerV12.txt
#runProc ${MYWDIR}/solA2v14 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v14.json MPSTEER=${MYWDIR}/mpSteerV12.txt
#runProc ${MYWDIR}/solA2v15 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v15.json MPSTEER=${MYWDIR}/mpSteerV12.txt
#runProc ${MYWDIR}/solA2v16 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v16.json MPSTEER=${MYWDIR}/mpSteerV12.txt
#runProc ${MYWDIR}/solA2v17 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v17.json MPSTEER=${MYWDIR}/mpSteerV12.txt
#runProc ${MYWDIR}/solA2v17N Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v17.json MPSTEER=${MYWDIR}/mpSteerV12.txt  # new type of global constraints
#runProc ${MYWDIR}/solA2v18 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v18.json MPSTEER=${MYWDIR}/mpSteerV12.txt  # same as v17 but w/o tangentialPolar_r0t1
#runProc ${MYWDIR}/solA2v18b0 Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v18.json MPSTEER=${MYWDIR}/mpSteerV12b0.txt  # same as v17 but w/o tangentialPolar_r0t1
#runProc ${MYWDIR}/solA2v18bpm Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v18.json MPSTEER=${MYWDIR}/mpSteerV12bpm.txt  # same as v17 but w/o tangentialPolar_r0t1

runProc ${MYWDIR}/tstDbg Config=/home/shahoian/alice/O2Align/doc/ITS_RB_v18.json MPSTEER=/data/devAlg/test1/sol0/mpSteer.txt  # 
