file(GLOB SOURCES "*.cpp")
add_library( HMC_Support SOURCES )
target_link_libraries( HMC_Support Actions Data)