file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_against_simple.f90 test_against_simple.f90)
target_link_libraries(test_against_simple.f90 vmec_lib utils_lib)
add_test(TestAgainstSIMPLE test_against_simple.f90)
set_tests_properties(TestAgainstSIMPLE PROPERTIES LABELS "quick")
