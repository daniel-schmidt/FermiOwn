add_executable( cmaketest "cmaketest.cpp" )

add_executable( Lattice_test "Lattice_test.cpp" )
target_link_libraries( Lattice_test Data )

add_executable( FieldScalar_test "FieldScalar_test.cpp" )
target_link_libraries( FieldScalar_test Data)

add_executable(Action_test "Action_test.cpp")
target_link_libraries( Action_test Actions Data)

add_executable(Metropolis_test "Metropolis_test.cpp")
target_link_libraries( Metropolis_test Actions Data UpdateAlgorithms)

add_executable( Integrator_test "Integrator_test.cpp" )
target_link_libraries( Integrator_test Actions Data HMC_Support )

add_executable(HMC_test "HMC_test.cpp")
target_link_libraries( HMC_test Actions Data UpdateAlgorithms)