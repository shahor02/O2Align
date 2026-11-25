#!/bin/bash

: ${MYWDIR:="/data0/devAlg/tests"}

ulimit -n 5000

: ${REMAP:=}
USECDB=""
echo "Imposed remapping for millepede.res -> alignment conversion: REMAP=$REMAP"

#export O2_DPL_MVBIAS=";mvbias.xyz[0]=0.033;mvbias.xyz[1]=-0.1825;mvbias.xyz[2]=0.0;mvbias.slopeX=-4.27e-05;mvbias.slopeY=-6.96e-05;"
#  copy initial geometry stuff 
if [[ -n ${STARTGEOM} ]] &&  [[ ${STARTGEOM} != "none" ]] ; then
    echo "Fetching initial geometry from ${STARTGEOM}"
    for det in GLO ITS TRD TOF MFT MCH ; do
	algName="${STARTGEOM}/${det}_alignment.root"
	if [[ -f "$algName" ]] ; then
	    cp ${algName} ./
            echo "Copied $algName to $(pwd)/"
	fi
	if [[ -d ${STARTGEOM}/${det} ]] ; then
	    echo "Sync ${STARTGEOM}/${det} to $(pwd)/"
	    rsync -av ${STARTGEOM}/${det} ./
	    [[ "0$USECDB" != "0" ]] && USECDB+=","
	    if [[ $det == "GLO" ]] ; then
		USECDB+="GLO/Config/GeometryAligned"
	    else 
		USECDB+="${det}/Calib/Align"
	    fi
	fi
    done 
fi
[[ -z $REMAP ]] && [[ "0$USECDB" != "0" ]] && REMAP="--condition-remap file://$(pwd)=$USECDB"
export REMAP=${REMAP:=}
echo "Final remapping for millepede.res -> alignment conversion: REMAP=$REMAP"


if [[ -z $PEDEDONE ]] ; then
   [[ -z $SKIPMILLESTEER ]] && runOPTL="-b --run"  Config=${Config:-} DataDir="N/A" OutType="MilleSteer" ${MYWDIR}/proc0.sh
   pede ${MPSTEER:=${MYWDIR}/mpSteer.txt} | tee pede.log    
fi

if [[ ! -f millepede.res ]] ; then
   echo "millepede.res file is not in the `pwd`"
   exit 1
fi

$O2Align_ROOT/share/macro/label_mille_results.py
OutType="MilleRes"
[[ -z $NOWEEKMODES ]] && OutType+=",WeakModesDiag"
runOPTL="-b --run"  Config=${Config:-} DataDir="N/A" OutType=$OutType ${MYWDIR}/proc0.sh
[[ -z $NOWEEKMODES ]] && root -b -q $O2Align_ROOT/share/macro/drawWeakModes.C+
    
TIMESTAMP="1746914681240"  $O2Align_ROOT/share/macro/makeCDB.sh

if [[ -z $SKIPCHECKRESID ]] ; then
  REMAP= ${MYWDIR}/batResRS.sh
  #EXTRARESOPT="--no-tree" ${MYWDIR}/batResRS.sh
  ${MYWDIR}/reResidAll.sh
fi
