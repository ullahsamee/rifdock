#!/bin/bash

# mkdir -p build_debug
cd build_debug
CXX=g++-6 CC=gcc-6 CMAKE_ROSETTA_PATH=/Users/brian/Documents/baker/rifdock/main CMAKE_FINAL_ROSETTA_PATH=/Users/brian/Documents/baker/rifdock/main/source/cmake/build_cxx11_omp_hdf5_debug cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j3 rif_dock_test rifgen
