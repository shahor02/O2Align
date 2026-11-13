#!/bin/bash


[[ -z $DESTCDB ]] && DESTCDB="."
[[ -z $ALGDIR ]] && ALGDIR="."
DESTCDB=`realpath $DESTCDB`
ALGDIR=`realpath $ALGDIR`

echo "Source alignment dir: $ALGDIR"
echo "Destination CCDB dir: $DESTCDB"

REMAP=""

for det in ITS TRD TOF MFT MCH ; do
  algName="${ALGDIR}/${det}_alignment.root"
  if [[ -f "$algName" ]] ; then
     [[ ! -d ${DESTCDB}/$det/Calib/Align ]] && mkdir -p ${DESTCDB}/${det}/Calib/Align
     cp $algName ${DESTCDB}/${det}/Calib/Align/snapshot.root
     [[ "0$REMAP" != "0" ]] && REMAP+=","
     REMAP+="${det}/Calib/Align"
  fi
done

if [[ "0$REMAP" != "0" ]] ; then
  REMAP="--condition-remap file://${DESTCDB}=${REMAP}"
  echo "Remapping CCDB to $REMAP"
else
  echo "No <det>_alignment.root file was found in the $ALGDIR"
  exit 1  
fi

if [[ -z $TIMESTAMP ]] ; then
  TIMESTAMP=$(($(date +%s%N)/1000000))
  echo "Timestamp for CCDB was not provided, set to now = $TIMESTAMP"
else
  echo "Using provided CCDB Timestamp $TIMESTAMP"
fi

[[ -f o2sim_geometry-aligned.root ]] && rm o2sim_geometry-aligned.root
o2-create-aligned-geometry-workflow -b --run $REMAP --condition-timestamp $TIMESTAMP

[[ ! -f o2sim_geometry-aligned.root ]] && { echo "failed to create aligned geometry" ; exit 1; }
[[ ! -d ${DESTCDB}/GLO/Config/GeometryAligned ]] && mkdir -p ${DESTCDB}/GLO/Config/GeometryAligned
cp o2sim_geometry-aligned.root ${DESTCDB}/GLO/Config/GeometryAligned/snapshot.root

[[ ! -f its_GeometryTGeo.root ]] && { echo "failed to create its_GeometryTGeo.root" ; exit 1; }
[[ ! -d ${DESTCDB}/ITS/Config/Geometry ]] && mkdir -p ${DESTCDB}/ITS/Config/Geometry
cp its_GeometryTGeo.root  ${DESTCDB}/ITS/Config/Geometry/snapshot.root

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
root -b -q $SCRIPT_DIR/algDump.C+\(\"o2sim_geometry-aligned.root\"\)


