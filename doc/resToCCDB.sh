#!/bin/bash
# Turn the millepede.res of the current directory into an ITS alignment with the local CCDB layout used by checkResid:
#   MilleRes (result.json, ITS_alignment.root)  ->  ITS/Calib/Align/snapshot.root
#   o2-create-aligned-geometry-workflow         ->  o2sim_geometry-aligned.root (+ GLO/Config/GeometryAligned/snapshot.root),
#                                                   its_GeometryTGeo.root (+ ITS/Config/Geometry/snapshot.root)
# The DOF config and the mean-vertex time slots must be those of the pede run that produced (or that the edit of) millepede.res
# is based on, otherwise the labels do not match.
# Usage (in the directory holding millepede.res, O2PDPSuite loaded, e.g. `o2d`):
#   Config=~/alice/O2Align/doc/ITS_RB_v19.json ~/alice/O2Align/doc/resToCCDB.sh
# Then checkResid can be run with CCDB="--condition-remap file://$(pwd)=GLO/Config/GeometryAligned".

: ${Config:?"set Config to the DOF config json of the solution"}
: ${MVJSON:="/data/devAlg/MV_B0BPBM.json"}  # AlignParams.MVTimeSlotsJson of the solution
: ${TimeStamp:=1746914681240}               # condition timestamp of the data (the default of proc0.sh for MilleRes)
O2ALIGN=${O2ALIGN:-~/alice/O2Align}

[[ -f millepede.res ]] || { echo "no millepede.res in $(pwd)"; exit 1; }

# O2 ships its own libO2Align.so: the local one must come first
export LD_LIBRARY_PATH=$O2ALIGN/install/lib:$LD_LIBRARY_PATH
export PATH=$O2ALIGN/install/bin:$PATH

runOPTL="-b --run" Config=$Config DataDir="N/A" OutType="MilleRes" TimeStamp=$TimeStamp \
  ConfString="AlignParams.dofConfigJson=$Config;AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=6;AlignParams.MVTimeSlotsJson=$MVJSON;" \
  $O2ALIGN/doc/proc0.sh < /dev/null > millres.log 2>&1
[[ -f ITS_alignment.root ]] || { echo "MilleRes failed, see millres.log"; exit 1; }

mkdir -p ITS/Calib/Align GLO/Config/GeometryAligned ITS/Config/Geometry
cp ITS_alignment.root ITS/Calib/Align/snapshot.root
# without an explicit timestamp the workflow queries the CCDB at the current time
o2-create-aligned-geometry-workflow --condition-remap file://$(pwd)=ITS/Calib/Align --condition-timestamp $TimeStamp -b --run < /dev/null > geo.log 2>&1
[[ -f o2sim_geometry-aligned.root ]] || { echo "aligned geometry failed, see geo.log"; exit 1; }
cp o2sim_geometry-aligned.root GLO/Config/GeometryAligned/snapshot.root
cp its_GeometryTGeo.root ITS/Config/Geometry/snapshot.root
echo "done: $(pwd)"
