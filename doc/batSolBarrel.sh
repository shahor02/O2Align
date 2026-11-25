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
  env "$@" $MYWDIR/sol.sh >& procSol.log
  popd
}

export STARTGEOM=${STARTGEOM:-"${MYWDIR}/best/IBEdge/"}
export DetList=${DetList:-"ITS-TRD,ITS-TOF,ITS-TRD-TOF"}

export MVSLOTS=${MVSLOTS:-"$MYWDIR/../data/MV_B0BPBM.json"}

#FIXME="AlignParams.groupModesRelToRoots=false;"  # needed for ITSEnvTRDTOFch3_RB.json but not for ch4

#export EXTRACONF="AlignParams.MVTimeSlotsJson=${MYWDIR}/../data/MV_B0BPBM.json;AlignParams.statInFile=${MYWDIR}/Mille_IBEdgeS4R/procStartColl.root;"
export EXTRACONF="$FIXME;AlignParams.statInFile=${MYWDIR}/Mille_IBEdgeS4R/procStatColl.root;AlignParams.corrType=1;"

#if [[ -n ${STARTGEOM} ]] &&  [[ ${STARTGEOM} != "none" ]] ; then
#  EXTRACONF+=";AlignParams.startFromMilleRes=${STARTGEOM}/millepede_labeled.res;"
#fi

export SKIPCHECKRESID=1


#runProc ${MYWDIR}/solB0A1 Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOF_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelV0.txt
#runProc ${MYWDIR}/solB1A1 Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOFch_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelV0.txt
#runProc ${MYWDIR}/solB1A2 Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOFch2_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelV0.txt
#runProc ${MYWDIR}/solB1A3 Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOFch3_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelHPLUT.txt
#runProc ${MYWDIR}/solB1A4 Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOFch4_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelFix.txt
runProc ${MYWDIR}/solB2A4   Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOFB0only_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelFixnoBM.txt
runProc ${MYWDIR}/solB2v4tR2 STARTMILLE=/data0/devAlg/tests/sag/v4tR2/millepede.res Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOFch4_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelFixnoBM.txt

