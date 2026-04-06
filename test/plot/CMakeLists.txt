file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../integration/input/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/input)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../integration/vmec/input/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/input)

add_executable(test_plot.x test_plot.f90)
target_link_libraries(test_plot.x myplot_lib utils_lib)
add_test(TestPlot test_plot.x)
set_tests_properties(TestPlot PROPERTIES LABELS "plot")

add_executable(plot_anti_sigma.x plot_anti_sigma.f90)
target_link_libraries(plot_anti_sigma.x myplot_lib
                                             rabe_lib
                                             helpers_lib)
add_test(PlotAntiSigma plot_anti_sigma.x)
set_tests_properties(PlotAntiSigma PROPERTIES LABELS "plot")

add_executable(plot_deviation_helical_anti_sigma.x plot_deviation_helical_anti_sigma.f90)
target_link_libraries(plot_deviation_helical_anti_sigma.x myplot_lib
                                                          rabe_lib
                                                          helpers_lib)
add_test(PlotDeviationHelicalAntiSigma plot_deviation_helical_anti_sigma.x)
set_tests_properties(PlotDeviationHelicalAntiSigma PROPERTIES LABELS "plot")

add_executable(plot_helical_anti_sigma_small_aspect.x plot_helical_anti_sigma_small_aspect.f90)
target_link_libraries(plot_helical_anti_sigma_small_aspect.x myplot_lib
                                                             rabe_lib
                                                             helpers_lib)
add_test(PlotHelicalAntiSigmaSmallAspect plot_helical_anti_sigma_small_aspect.x)
set_tests_properties(PlotHelicalAntiSigmaSmallAspect PROPERTIES LABELS "plot")

add_executable(plot_deviation_poloidal_anti_sigma.x plot_deviation_poloidal_anti_sigma.f90)
target_link_libraries(plot_deviation_poloidal_anti_sigma.x myplot_lib
                                                          rabe_lib
                                                          helpers_lib)
add_test(PlotDeviationPoloidalAntiSigma plot_deviation_poloidal_anti_sigma.x)
set_tests_properties(PlotDeviationPoloidalAntiSigma PROPERTIES LABELS "plot")

add_executable(plot_local_drift.x plot_local_drift.f90)
target_link_libraries(plot_local_drift.x myplot_lib
                                             rabe_lib
                                             helpers_lib)
add_test(PlotLocalDrift plot_local_drift.x)
set_tests_properties(PlotLocalDrift PROPERTIES LABELS "plot")

add_executable(plot_nabla_s_compare.x plot_nabla_s_compare.f90)
target_link_libraries(plot_nabla_s_compare.x myplot_lib
                                              rabe_lib
                                              helpers_lib)
add_test(PlotNablaSCompare plot_nabla_s_compare.x)
set_tests_properties(PlotNablaSCompare PROPERTIES LABELS "current")
