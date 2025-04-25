file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_against_neo2_field.x test_against_neo2_field.f90)
target_link_libraries(test_against_neo2_field.x rabe_lib)
add_test(TestAgainstNeo2Field test_against_neo2_field.x)
set_tests_properties(TestAgainstNeo2Field PROPERTIES LABELS "slow")

add_executable(test_fieldline.x test_fieldline.f90)
target_link_libraries(test_fieldline.x rabe_lib)
add_test(TestFieldline test_fieldline.x)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/helpers)
add_subdirectory(helpers)

add_executable(test_integration_along_fieldline.x test_integration_along_fieldline.f90)
target_link_libraries(test_integration_along_fieldline.x rabe_lib helpers_lib)
add_test(TestIntegrationFieldline test_integration_along_fieldline.x)
