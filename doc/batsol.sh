#!/bin/bash


mkdir solA2v5bow 
pushd solA2v5bow
/data0/devAlg/tests/batResRS.sh
popd

pushd /data0/devAlg/tests/sol1/savITSBarFix
/data0/devAlg/tests/batResRS.sh
popd

