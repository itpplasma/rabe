file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_against_neo2.x test_against_neo2.f90)
target_link_libraries(test_against_neo2.x rabe_lib)
add_test(TestAgainstNeo2 test_against_neo2.x)
set_tests_properties(TestAgainstNeo2 PROPERTIES LABELS "slow")

add_executable(test_fieldline.x test_fieldline.f90)
target_link_libraries(test_fieldline.x rabe_lib)
add_test(TestFieldline test_fieldline.x)

add_library(mock_field STATIC
  mock_field.f90
)
target_link_libraries(mock_field PRIVATE rabe_lib)

add_executable(test_integration_along_fieldline.x test_integration_along_fieldline.f90)
target_link_libraries(test_integration_along_fieldline.x rabe_lib mock_field)
add_test(TestIntegrationFieldline test_integration_along_fieldline.x)
