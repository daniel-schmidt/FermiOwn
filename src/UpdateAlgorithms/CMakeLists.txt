add_subdirectory("HMC_Support")
file(GLOB SOURCES "*.cpp")
add_library( UpdateAlgorithms ${SOURCES} )
target_link_libraries( UpdateAlgorithms Actions Data)