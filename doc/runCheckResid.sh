#!/bin/bash
export DPL_REPORT_PROCESSING=1

: ${PVCONSTR:="true"}

[[ -n $PVSTAT ]] && [[ $PVSTAT == "0" ]] && PVCONSTR="false"

echo "CCDB = $CCDB, PIPELINES = $PIPELINES CONFIG=$CONFIG NTHREADS=$NTHREADS PVSTAT=$PVSTAT PVCONSTR=$PVCONSTR"

GLOSET="--shm-segment-size 16000000000  --timeframes-rate-limit ${PIPELINES:-4} --timeframes-rate-limit-ipcid $(($(date +%s%3N) % 1234568))"
PVCONF=";pvertexer.useMeanVertexConstraint=false;pvertexer.meanVertexExtraErrSelection=0.2;pvertexer.iniScale2=100;pvertexer.acceptableScale2=10.;"

cmd="o2-reader-driver-workflow --max-tf ${MAXTF:--1}  $GLOSET --hbfutils-config o2_tfidinfo.root,upstream | \
o2-check-resid-workfow $CCDB $GLOSET $EXTRARESOPT --nthreads ${NTHREADS:-10} --hbfutils-config o2_tfidinfo.root,upstream --pipeline check-resid:${PIPELINES:-1} \
--configKeyValues \"checkresid.addPVAsCluster=${PVCONSTR};${CONFIG};${PVCONF};$EXTRAOPT\"  --run -b"

echo "Running $cmd"
eval $cmd

ls -rlat *.root


if [[ -f checkResid_t0_hman.root ]] ; then
  echo "Merging $(ls checkResid_t*_hman.root)"
  hadd -f checkResid_hman.root checkResid_t*_hman.root
  rm checkResid_t*_hman.root
fi

if [[ -f checkResid_t0.root ]] ; then
  echo "Merging $(ls checkResid_t*.root)"
  hadd -f checkResid.root checkResid_t*.root
  rm checkResid_t*.root
fi
