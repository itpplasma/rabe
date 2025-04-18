file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_against_neo2.x test_against_neo2.f90)
target_link_libraries(test_against_neo2.x rabe_lib)
add_test(TestAgainstNeo2 test_against_neo2.x)
