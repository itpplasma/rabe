file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_against_simple.x test_against_simple.f90)
target_link_libraries(test_against_simple.x vmec_lib utils_lib)
add_test(TestAgainstSIMPLE test_against_simple.x)
set_tests_properties(TestAgainstSIMPLE PROPERTIES LABELS "quick")

add_executable(test_against_neo_qa.x test_against_neo_qa.f90)
target_link_libraries(test_against_neo_qa.x vmec_lib utils_lib)
add_test(TestAgainstNEO_QA test_against_neo_qa.x)
set_tests_properties(TestAgainstNEO_QA PROPERTIES LABELS "quick")

add_executable(test_against_neo_qh.x test_against_neo_qh.f90)
target_link_libraries(test_against_neo_qh.x vmec_lib utils_lib)
add_test(TestAgainstNEO_QH test_against_neo_qh.x)
set_tests_properties(TestAgainstNEO_QH PROPERTIES LABELS "quick")

add_executable(test_surface_area.x test_surface_area.f90)
target_link_libraries(test_surface_area.x vmec_lib utils_lib)
add_test(TestSurfaceArea test_surface_area.x)
set_tests_properties(TestSurfaceArea PROPERTIES LABELS "quick")
