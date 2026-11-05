#!/bin/bash

: ${PVSTAT:="1"}

: ${OUTNAME:="checkResid"}

: ${DTB0:="checkResid_LHC25ab_562862"}
: ${LBB0:="pp B0(LHC25ab/562862)"}

: ${DTBP:="checkResid_LHC25ae"}
: ${LBBP:="OO B+(LHC25ae/00564356)"}

echo "PVSTAT=$PVSTAT"
echo "OutName=$OUTNAME"
echo "DataB0=$DTB0 LabelB0=$LBB0"
echo "DataBP=$DTBP LabelBP=$LBBP"

cmd="o2-check-resid-workfow --draw-external-only --postproc-external-only --configKeyValues \
\"checkresid.outname=${OUTNAME}_PV${PVSTAT};\
checkresid.savePostprocessed=true;\
checkresid.ext_hm_list=\
/data0/devAlg/tests/v10n/${DTBP}_hman_AlgV10n_PV${PVSTAT}.root,/data0/devAlg/tests/v10n/${DTB0}_hman_AlgV10n_PV${PVSTAT}.root,\
${DTBP}_hman_AlgRS_PV${PVSTAT}.root,${DTB0}_hman_AlgRS_PV${PVSTAT}.root\
;\
checkresid.ext_leg_list=\
V10n ${LBBP},V10n ${LBB0},\
MP2 ${LBBP},MP2 ${LBB0}\
;\
checkresid.ext_color_list=\
kGreen+2,kBlack,\
kBlue,kRed;
checkresid.ext_marker_list=\
20,22,\
24,26;\
\" -b --run
"

echo "Running $cmd"
eval $cmd
