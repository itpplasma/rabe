add_executable(test_find_extrema.x test_find_extrema.f90)
target_link_libraries(test_find_extrema.x rabe_lib)
add_test(TestFindExtrema test_find_extrema.x)
set_tests_properties(TestFindExtrema PROPERTIES LABELS "quick")

add_executable(test_integrate_1d.x test_integrate_1d.f90)
target_link_libraries(test_integrate_1d.x rabe_lib)
add_test(TestIntegrate1d test_integrate_1d.x)
set_tests_properties(TestIntegrate1d PROPERTIES LABELS "quick")

add_executable(test_integrate_1d_substituted.x test_integrate_1d_substituted.f90)
target_link_libraries(test_integrate_1d_substituted.x rabe_lib)
add_test(TestIntegrate1dSubstituted test_integrate_1d_substituted.x)
set_tests_properties(TestIntegrate1dSubstituted PROPERTIES LABELS "quick")

add_executable(test_integrate_wrapper.x test_integrate_wrapper.f90)
target_link_libraries(test_integrate_wrapper.x rabe_lib)
add_test(TestIntegrateWrapper test_integrate_wrapper.x)
set_tests_properties(TestIntegrateWrapper PROPERTIES LABELS "quick")

add_executable(test_real_ft.x test_real_ft.f90)
target_link_libraries(test_real_ft.x rabe_lib)
add_test(TestRealFt test_real_ft.x)
set_tests_properties(TestRealFt PROPERTIES LABELS "quick")

add_executable(test_misc.x test_misc.f90)
target_link_libraries(test_misc.x rabe_lib)
add_test(TestMisc test_misc.x)
set_tests_properties(TestMisc PROPERTIES LABELS "quick")

add_executable(plot_misc.x plot_misc.f90)
target_link_libraries(plot_misc.x rabe_lib
                                  myplot_lib)
add_test(PlotMisc plot_misc.x)
set_tests_properties(PlotMisc PROPERTIES LABELS "plot")
