add_executable(test_find_extrema.x test_find_extrema.f90)
target_link_libraries(test_find_extrema.x rabe_lib)

add_test(TestFindExtrema test_find_extrema.x)