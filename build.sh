#!/bin/sh
module purge
module load cmake gcc/8.3.1 mkl/2019.5 hpcx/2.13.0

USE_SCOREP=0
REBUILD=0

if [[ $REBUILD -eq 1 ]]; then
    rm -rf build
    mkdir -p build
fi
cd build

if [[ $USE_SCOREP -eq 1 ]]; then
    module load scorep/8.0-hpcx-2.13.0
    if [[ $REBUILD -eq 1 ]]; then
        ecbuild_ROOT=$HOME/applications/modules/ecbuild/3.7.1/ fiat_ROOT=$HOME/applications/modules/fiat/1.1.0/ CC=scorep-mpicc FC=scorep-mpif90 cmake -DCMAKE_INSTALL_PREFIX=$HOME/applications/modules/ectrans/1.2.0-scorep -DENABLE_OMP:BOOL=OFF .. 2>&1 | tee configure.log
    fi
else
    if [[ $REBUILD -eq 1 ]]; then
        ecbuild_ROOT=$HOME/applications/modules/ecbuild/3.7.1/ fiat_ROOT=$HOME/applications/modules/fiat/1.1.0/ CC=mpicc FC=mpif90 cmake -DCMAKE_INSTALL_PREFIX=$HOME/applications/modules/ectrans/1.2.0 .. 2>&1 | tee configure.log
    fi
fi

make install 2>&1 | tee make.log
