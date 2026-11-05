#!/bin/bash

#export O2_DPL_MVBIAS=";mvbias.xyz[0]=0.033;mvbias.xyz[1]=-0.1825;mvbias.xyz[2]=0.0;mvbias.slopeX=-4.27e-05;mvbias.slopeY=-6.96e-05;"

if [[ -z $PEDEDONE ]] ; then
    runOPTL="-b --run"  Config=${Config:-} DataDir="N/A" OutType="MilleSteer" /data0/devAlg/tests/proc0.sh
    pede ${MPSTEER:=/data0/devAlg/tests/mpSteer.txt} | tee pede.log
fi
$O2Align_ROOT/share/macro/label_mille_results.py
runOPTL="-b --run"  Config=${Config:-} DataDir="N/A" OutType="MilleRes,WeakModesDiag" /data0/devAlg/tests/proc0.sh
root -b -q ~/alice/O2Align/macro/drawWeakModes.C+
TIMESTAMP="1746914681240"  ~/alice/mymacrosO2/makeCDB.sh

/data0/devAlg/tests/batResRS.sh
#EXTRARESOPT="--no-tree" /data0/devAlg/tests/batResRS.sh
/data0/devAlg/tests/reResidAll.sh

