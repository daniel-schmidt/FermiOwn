cmake_minimum_required( VERSION 2.8 )
project ( FermiOwn )
include(CTest)

set( FermiOwn_VERSION_MAJOR 0 )
set( FermiOwn_VERSION_MINOR 1 )

set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11 -lm -Wall")
#set( EIGEN3_INCLUDE_DIR "/work/libs/eigen3.2.7/install/include/eigen3" )

include_directories($ENV{EIGEN3_INCLUDE_DIR})
include_directories("src/Math")
include_directories("src/Data")
include_directories("src/Actions")
include_directories("src/UpdateAlgorithms")
include_directories("src/UpdateAlgorithms/HMC_Support")

add_subdirectory("src")
add_subdirectory("test")
add_subdirectory("apps")