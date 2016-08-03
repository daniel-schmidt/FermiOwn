add_executable( EigenBool_test "EigenBool_test.cpp" )

add_executable( FieldBoolean_test "FieldBoolean_test.cpp" )
target_link_libraries( FieldBoolean_test Data)

# testing the Slac Operator Matrix implementation in 1D for sizes 2 to 5
# add_executable(SlacOperatorMatrix_test "SlacOperatorMatrix_test.cpp")
# target_link_libraries( SlacOperatorMatrix_test Actions)

#add_test( SlacOperatorMatrix_testRuns SlacOperatorMatrix_test 3 )

#macro (do_slac_test arg result)
#  add_test (SlacOperatorMatrix_test${arg} SlacOperatorMatrix_test ${arg})
#  set_tests_properties ( SlacOperatorMatrix_test${arg}
#    PROPERTIES PASS_REGULAR_EXPRESSION ${result})
#endmacro (do_slac_test)

#do_slac_test( 2, "Size 2 works!" )
#do_slac_test( 3, "Size 3 works!" )
#do_slac_test( 4, "Size 4 works!" )
#do_slac_test( 5, "Size 5 works!" )

# TODO: refine regex to find whole det, not just decimal places...
#add_test( SlacOperatorMatrix_testDet SlacOperatorMatrix_test 4 )
#set_tests_properties( SlacOperatorMatrix_testDet PROPERTIES PASS_REGULAR_EXPRESSION "69524" )

add_executable( WoodburyMatrix_test "WoodburyMatrix_test.cpp" )
target_link_libraries( WoodburyMatrix_test Math )

add_executable( ConfigGenerator_test "ConfigGenerator_test.cpp" )
target_link_libraries( ConfigGenerator_test Actions )

add_executable( WeightFunction_test "WeightFunction_test.cpp" )
target_link_libraries( WeightFunction_test Actions )

add_executable( DSlashUpdater_test "DSlashUpdater_test.cpp" )
target_link_libraries( DSlashUpdater_test Actions )