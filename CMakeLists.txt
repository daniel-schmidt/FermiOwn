cmake_minimum_required( VERSION 2.8 )
project ( FermiOwn )
include(CTest)

set( FermiOwn_VERSION_MAJOR 0 )
set( FermiOwn_VERSION_MINOR 1 )

set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11 -lm -Wall -pg -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native -std=c++11 -lm -Wall -O0 -g")

include_directories($ENV{EIGEN3_INCLUDE_DIR})
include_directories("src/Math")
include_directories("src/Data")
include_directories("src/Actions")
include_directories("src/UpdateAlgorithms")
include_directories("FermiOwnCommon/src/Actions")
include_directories("FermiOwnCommon/src/Algorithms")
include_directories("FermiOwnCommon/src/Data")
include_directories("FermiOwnCommon/src/Math")
include_directories("FermiOwnCommon/src/System")


add_subdirectory("FermiOwnCommon")
add_subdirectory("src")
add_subdirectory("test")
add_subdirectory("apps")