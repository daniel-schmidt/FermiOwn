Simulation of the irreducible Thirring model after Fierz transformation, rewritten in a fermion bag like way.

Tests are provided if sensible...

shared functionality is included in the submodule FermiOwnCommon, clone with git clone --recursive to include it, or after clone:
git submodule init
git submodule update --recursive

in console: 

export CPLUS_INCLUDE_PATH="/work/libs/eigen3.2.7/install/include/eigen3"
mkdir build
cd build
cmake ../
make

switch build types with cmake -DCMAKE_BUILD_TYPE=Debug or -DCMAKE_BUILD_TYPE=Release

profiling:

* gprof:
apps/FermiOwnGraph 2 3 200 1000 100 0 1 > ref-1.log 2> ref-1.dat
gprof apps/FermiOwnGraph | gprof2dot -sw --skew 0.01 | dot -Tpng -o output.png

* operf:
operf --callgraph apps/FermiOwnGraph 2 3 200 1000 100 0 1
opreport -cgf | gprof2dot -f oprofile -sw --skew 0.01 | dot -Tpng -o out.png