add_executable( CliffordAlgebra_test "CliffordAlgebra_test.cpp" )
add_test( CliffordAlgebra_testAlgebra CliffordAlgebra_test )
set_tests_properties( CliffordAlgebra_testAlgebra PROPERTIES PASS_REGULAR_EXPRESSION  "3 matrices squared to unity." )
set_tests_properties( CliffordAlgebra_testAlgebra PROPERTIES PASS_REGULAR_EXPRESSION  "3  correct products." )
set_tests_properties( CliffordAlgebra_testAlgebra PROPERTIES PASS_REGULAR_EXPRESSION  "3 commutator relations worked." )

add_executable( EigenBool_test "EigenBool_test.cpp" )

add_executable( Lattice_test "Lattice_test.cpp" )
target_link_libraries( Lattice_test Data )

add_executable( Field_test "Field_test.cpp" )
target_link_libraries( Field_test Data)
add_test( Field_testRuns Field_test )
set_tests_properties (Field_testRuns PROPERTIES PASS_REGULAR_EXPRESSION "Initialization with 0 and 1 works!")
set_tests_properties (Field_testRuns PROPERTIES PASS_REGULAR_EXPRESSION "Initialization by copy works!")
set_tests_properties (Field_testRuns PROPERTIES PASS_REGULAR_EXPRESSION "Copy assignment works!")

add_executable(Action_test "Action_test.cpp")
target_link_libraries( Action_test Actions Data)

add_executable(Metropolis_test "Metropolis_test.cpp")
target_link_libraries( Metropolis_test Actions Data UpdateAlgorithms)

add_executable( Integrator_test "Integrator_test.cpp" )
target_link_libraries( Integrator_test Actions Data HMC_Support )

add_executable(HMC_test "HMC_test.cpp")
target_link_libraries( HMC_test Actions Data UpdateAlgorithms HMC_Support)

# testing the Slac Operator Matrix implementation in 1D for sizes 2 to 5
add_executable(SlacOperatorMatrix_test "SlacOperatorMatrix_test.cpp")
target_link_libraries( SlacOperatorMatrix_test Actions)

add_test( SlacOperatorMatrix_testRuns SlacOperatorMatrix_test 3 )

macro (do_slac_test arg result)
  add_test (SlacOperatorMatrix_test${arg} SlacOperatorMatrix_test ${arg})
  set_tests_properties ( SlacOperatorMatrix_test${arg}
    PROPERTIES PASS_REGULAR_EXPRESSION ${result})
endmacro (do_slac_test)

do_slac_test( 2, "Size 2 works!" )
do_slac_test( 3, "Size 3 works!" )
do_slac_test( 4, "Size 4 works!" )
do_slac_test( 5, "Size 5 works!" )

# TODO: refine regex to find whole det, not just decimal places...
add_test( SlacOperatorMatrix_testDet SlacOperatorMatrix_test 4 )
set_tests_properties( SlacOperatorMatrix_testDet PROPERTIES PASS_REGULAR_EXPRESSION "69524" )