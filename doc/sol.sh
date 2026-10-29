#!/bin/bash

MYWDIR="/data/devAlg"

#export O2_DPL_MVBIAS=";mvbias.xyz[0]=0.033;mvbias.xyz[1]=-0.1825;mvbias.xyz[2]=0.0;mvbias.slopeX=-4.27e-05;mvbias.slopeY=-6.96e-05;"

if [[ -z $PEDEDONE ]] ; then
  runOPTL="-b --run"  Config=${Config:-} DataDir="N/A" OutType="MilleSteer" ${MYWDIR}/proc0.sh
  pede ${MPSTEER:=$MYWDIR/mpSteer.txt} | tee pede.log
fi

$O2Align_ROOT/share/macro/label_mille_results.py
runOPTL="-b --run"  Config=${Config:-} DataDir="N/A" OutType="MilleRes,WeakModesDiag" ${MYWDIR}/proc0.sh
root -b -q ~/alice/O2Align/macro/drawWeakModes.C+
TIMESTAMP="1746914681240"  ~/alice/mymacrosO2/makeCDB.sh

CCDB="--condition-remap file://$(pwd)=GLO/Config/GeometryAligned" ${MYWDIR}/batResRS.sh
#EXTRARESOPT="--no-tree" ${MYWDIR}/batResRS.sh
${MYWDIR}/reResidAll.sh

