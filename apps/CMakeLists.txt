add_executable(FermiOwnGraph "FermiOwnGraph.cpp")
target_link_libraries( FermiOwnGraph UpdateAlgorithms SystemCommon )

add_executable(ExactPartitionSum "ExactPartitionSum.cpp")
target_link_libraries( ExactPartitionSum Actions )