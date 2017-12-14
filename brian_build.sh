#!/bin/bash

#rm -r build/*
cd build
CXX=g++-6 CC=gcc-6 CMAKE_ROSETTA_PATH=/Users/brian/Documents/baker/rifdock/main CMAKE_FINAL_ROSETTA_PATH=/Users/brian/Documents/baker/rifdock/main/source/cmake/build_cxx11_omp_hdf5 cmake .. -DCMAKE_BUILD_TYPE=Release
make -j1 rif_dock_test rifgen 
