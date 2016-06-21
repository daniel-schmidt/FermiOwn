add_executable(FermiOwnGraph "FermiOwnGraph.cpp")
target_link_libraries( FermiOwnGraph Actions Data UpdateAlgorithms )

add_executable(ExactPartitionSum "ExactPartitionSum.cpp")
target_link_libraries( ExactPartitionSum Actions )