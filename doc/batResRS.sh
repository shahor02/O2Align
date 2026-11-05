#!/bin/bash

procDir()  # args: directory name
{
  if [[ $# -lt 2 ]]; then
      echo "procDir expects arguments: <dataDir> <jobName>"
    exit 1
  fi
  ddir=$1
  DNAME=$2
  pushd $1
  echo "Processing in dir $ddir Name=$DNAME PVSTAT=$PVSTAT"
  EXTRARESOPT=$EXTRARESOPT CONFIG="$CONFIG;checkresid.doIBOBresiduals=true;" CCDB=$REMAP PIPELINES=8 NTHREADS=4 MAXTF=${MAXTF} PVSTAT=${PVSTAT}  $script
  if [[ -f checkResid.root ]] ; then
      cmd="mv checkResid.root ${DEST}/checkResid_${DNAME}_${PREF}_PV${PVSTAT}.root"
      echo "running: $cmd"
      eval $cmd
  else
      echo "Output residuals tree for $DNAME is absent in $ddir"
  fi
  if [[ -f checkResid_hman.root ]] ; then
      cmd="mv checkResid_hman.root ${DEST}/checkResid_${DNAME}_hman_${PREF}_PV${PVSTAT}.root"
      echo "running: $cmd"
      eval $cmd
  else
      echo "Output residuals hman for $DNAME is absent in $ddir"
  fi
  popd
}

#EXTRARESOPT="--no-tree --draw-report"
script="/data0/devAlg/tests/runCheckResid.sh"

: ${PVSTAT:="1"}

: ${CONFIG:="checkresid.nBinsPhi=100;checkresid.nBinsZ=100;checkresid.nBinsTgl=60;"}
#: ${CONFIG:="checkresid.nBinsPhi=100;checkresid.nBinsZ=100;checkresid.nBinsTgl=60;checkresid.skipLr=1;"}

: ${DEST:=`pwd`}

: ${REMAP:="--condition-remap file://$(pwd)=GLO/Config/GeometryAligned"}

: ${PREF:="AlgRS"}
#REMAP='--condition-remap file:///data0/ALG25/Pede_v10n=GLO/Config/GeometryAligned'

MAXTF=4000
procDir /data0/devAlg/data/B0 LHC25ab_562862
MAXTF=12000
procDir /data0/devAlg/data/BPlus LHC25ab_562837
procDir /data0/devAlg/data/BMinus LHC25ah_564644
#procDir /dataHD/OO/LHC25ae_run00564356 LHC25ae
