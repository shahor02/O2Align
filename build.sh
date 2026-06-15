#!/bin/bash

# to be run from the Alig dir.
[[ -z $O2_ROOT ]] && alienv load O2PDPSuite/latest Millepede-II/latest ; source ~/bin/o2d.sh ;
export ALGO2_ROOT=`pwd`

[[ ! -d $ALGO2_ROOT/build ]] && mkdir $ALGO2_ROOT/build
[[ ! -d $ALGO2_ROOT/install ]] && mkdir $ALGO2_ROOT/install

cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$ALGO2_ROOT/install  -G Ninja
cd $ALGO2_ROOT/build
ninja -j 20 install
