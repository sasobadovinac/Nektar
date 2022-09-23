#!/bin/bash -x

[[ $OS_VERSION != "osx" ]] && ccache -s && ccache -M 5G

if [[ $BUILD_TYPE == "default" ]]; then
    BUILD_OPTS="-DCMAKE_BUILD_TYPE=Release \
        -DNEKTAR_TEST_ALL=ON \
        -DNEKTAR_ERROR_ON_WARNINGS=OFF"
elif [[ $BUILD_TYPE == "full" ]]; then
    BUILD_OPTS="-DCMAKE_BUILD_TYPE:STRING=Debug \
        -DNEKTAR_FULL_DEBUG:BOOL=ON \
        -DNEKTAR_TEST_ALL:BOOL=ON \
        -DNEKTAR_USE_ARPACK:BOOL=ON \
        -DNEKTAR_USE_FFTW:BOOL=ON \
        -DNEKTAR_USE_MPI:BOOL=ON \
        -DNEKTAR_USE_SCOTCH:BOOL=ON \
        -DNEKTAR_USE_PETSC:BOOL=ON \
        -DNEKTAR_USE_HDF5:BOOL=ON \
        -DNEKTAR_USE_MESHGEN:BOOL=ON \
        -DNEKTAR_USE_CCM:BOOL=ON \
        -DNEKTAR_CCMIO_URL=https://www.nektar.info/ccmio/libccmio-2.6.1.tar.gz \
        -DNEKTAR_USE_CWIPI:BOOL=ON \
        -DNEKTAR_USE_VTK:BOOL=ON \
        -DNEKTAR_BUILD_PYTHON:BOOL=ON \
        -DNEKTAR_TEST_USE_HOSTFILE=ON \
        -DNEKTAR_ERROR_ON_WARNINGS=OFF"
    if [[ $BUILD_SIMD == "avx2" ]]; then
        BUILD_OPTS="$BUILD_OPTS -DNEKTAR_ENABLE_SIMD_AVX2:BOOL=ON"
    elif [[ $BUILD_SIMD == "avx512" ]]; then
        BUILD_OPTS="$BUILD_OPTS -DNEKTAR_ENABLE_SIMD_AVX512:BOOL=ON"
    fi
    if [[ $BUILD_SIMD == "avx2" ]]; then
        BUILD_OPTS="$BUILD_OPTS -DNEKTAR_ENABLE_SIMD_AVX2"
    fi
fi

# Custom compiler
if [[ $BUILD_CC != "" ]]; then
   BUILD_OPTS="$BUILD_OPTS -DCMAKE_C_COMPILER=${BUILD_CC}"
fi
if [[ $BUILD_CXX != "" ]]; then
   BUILD_OPTS="$BUILD_OPTS -DCMAKE_CXX_COMPILER=${BUILD_CXX}"
fi
if [[ $BUILD_FC != "" ]]; then
   BUILD_OPTS="$BUILD_OPTS -DCMAKE_Fortran_COMPILER=${BUILD_FC}"
fi

rm -rf build && mkdir -p build && (cd build && cmake -G 'Unix Makefiles' $BUILD_OPTS ..) && \
    make -C build -j $NUM_CPUS all 2>&1 && \
    make -C build -j $NUM_CPUS install && \
    (cd build && ctest -j $NUM_CPUS --output-on-failure)

exit_code=$?
if [[ $exit_code -ne 0 ]]
then
    [[ $OS_VERSION != "osx" ]] && rm -rf build/dist
    exit $exit_code
fi
