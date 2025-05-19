add_executable(test_find_extrema.x test_find_extrema.f90)
target_link_libraries(test_find_extrema.x rabe_lib)
add_test(TestFindExtrema test_find_extrema.x)

add_executable(test_integrate_1d.x test_integrate_1d.f90)
target_link_libraries(test_integrate_1d.x rabe_lib)
add_test(TestIntegrate1d test_integrate_1d.x)

add_executable(test_integrate_wrapper.x test_integrate_wrapper.f90)
target_link_libraries(test_integrate_wrapper.x rabe_lib)
add_test(TestIntegrateWrapper test_integrate_wrapper.x)
