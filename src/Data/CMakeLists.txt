file(GLOB SOURCES "*.cpp")
add_library( Data ${SOURCES} )
target_link_libraries( Data DataCommon )