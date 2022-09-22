#!/bin/bash

set -e

if [ "$#" -ne 3 ]; then
    echo "Usage: test.sh /path/to/nektar++/cmake/dir ncpus /usr/bin/cmake"
    echo ""
    echo "Given a path to the Nektar++ CMake files, this script compiles and runs a simulation"
    echo "using ncpus, using the CMake executable given in the third parameter."
    exit 1
fi

# Create a build directory and compile against Nektar++
rm -rf build && mkdir build && cd build
$3 -DNektar++_DIR=$1 ..
make install

# Run test case in parallel.
cd ..
export OMPI_MCA_btl_vader_single_copy_mechanism=none
if [[ $EUID -eq 0 ]]; then
    mpiargs="-n $2 --allow-run-as-root"
else
    mpiargs="-n $2"
fi

test_output=`mpirun $mpiargs ./build/dist/ExampleSolver sample-laplace.xml | grep "L 2 error" | awk '{print ($7 < 1e-10)}'`
if [ "$test_output" -eq 1 ]; then
    echo "Test passed tolerance"
    exit 0
fi

# Test failed tolerance
echo "Test failed tolerance"
exit 1
