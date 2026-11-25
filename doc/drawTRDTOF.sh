#!/bin/bash

: ${REF:="/data0/devAlg/tests/v10n"}
: ${CUR:="."}
: ${RTP:="unb"}

root -b -q ${O2Align_ROOT}/share/macro/trdTofResid.C+\(\"./coll_LutPt08ITSSEED/_data0_devAlg_data_B0/debug_alg.root\",\"resTRDTOF_${RTP}_b0_hman.root\",\"o2sim_geometry-aligned.root\",0\)
root -b -q ${O2Align_ROOT}/share/macro/trdTofResid.C+\(\"./coll_LutPt08ITSSEED/_data0_devAlg_data_BPlus/debug_alg.root\",\"resTRDTOF_${RTP}_bp_hman.root\",\"o2sim_geometry-aligned.root\",0\)
root -b -q ${O2Align_ROOT}/share/macro/trdTofResid.C+\(\"./coll_LutPt08ITSSEED/_data0_devAlg_data_BMinus/debug_alg.root\",\"resTRDTOF_${RTP}_bm_hman.root\",\"o2sim_geometry-aligned.root\",0\)


root -l -b -q ${O2Align_ROOT}/share/macro/trdTofResid.C+ -e drawTRDTOFResid\(\
\""checkresid.ext_hm_list=\
${REF}/resTRDTOF_${RTP}_b0_hman.root,${REF}/resTRDTOF_${RTP}_bp_hman.root,${REF}/resTRDTOF_${RTP}_bm_hman.root,\
${CUR}/resTRDTOF_${RTP}_b0_hman.root,${CUR}/resTRDTOF_${RTP}_bp_hman.root,${CUR}/resTRDTOF_${RTP}_bm_hman.root;\
checkresid.ext_leg_list=\
v10n B0,v10n B+,v10n B-,\
mp2b B0,mp2b B+,mp2b B-;\
checkresid.ext_color_list=\
kBlack,kGreen+2,kYellow+2,\
kRed,kBlue,kViolet;\
checkresid.ext_marker_list=\
20,21,22,\
24,25,26;"\",\
\""cmp_${RTP}.pdf"\"\)
