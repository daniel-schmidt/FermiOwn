file(GLOB SOURCES "*.cpp")
add_library( Actions ${SOURCES} )
target_link_libraries( Actions Data Math )