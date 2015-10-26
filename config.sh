#!/bin/bash

cmake \
    -DCMAKE_CXX_COMPILER=`which g++` \
    -DCMAKE_C_COMPILER=`which gcc` \
    -DVTK_DIR=/work/visit/warp-visit/visit-deps/visit/vtk/6.1.0/x86_64r/lib/cmake/vtk-6.1 \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS_RELEASE="-Wall -Wextra -DNDEBUG -Ofast -march=native -mavx -ffast-math -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-local-typedefs -fopt-info-optimized" \
    -DCMAKE_CXX_FLAGS_DEBUG="-g3 -Wall -Wextra" \
    ..

