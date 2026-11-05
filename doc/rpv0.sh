#!/bin/bash

export PVSTAT="1"
export DTBP="checkResid_LHC25ab_562837"
export LBBP="pp B+(LHC25ab/562837)"
export OUTNAME="checkResidX"

pushd /data0/devAlg/tests/solA2v18
/data0/devAlg/tests/batResRS.sh
/data0/devAlg/tests/reResidAll.sh
popd

exit
 
pushd /data0/devAlg/tests/v10n
./batResRS.sh
popd

pushd /data0/devAlg/tests/solA2v19
/data0/devAlg/tests/batResRS.sh
/data0/devAlg/tests/reResidAll.sh
popd

pushd /data0/devAlg/tests/solA2v19cure/cure012
/data0/devAlg/tests/batResRS.sh
/data0/devAlg/tests/reResidAll.sh
popd


