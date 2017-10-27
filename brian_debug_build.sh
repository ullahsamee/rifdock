#!/bin/bash

mkdir -p build_debug
cd build_debug
CMAKE_ROSETTA_PATH=/home/bcov/rifdock/rosetta_for_rif/master/ cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j20 rif_dock_test rifgen
