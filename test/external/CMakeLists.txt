file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(plot_deviation_landreman_paul_helical.x plot_deviation_landreman_paul_helical.f90)
target_link_libraries(plot_deviation_landreman_paul_helical.x myplot_lib
                                                          rabe_lib
                                                          helpers_lib)
add_test(PlotDeviationLandremanPaulHelical plot_deviation_landreman_paul_helical.x)
set_tests_properties(PlotDeviationLandremanPaulHelical PROPERTIES LABELS "external")

add_executable(plot_deviation_landreman_paul_axi.x plot_deviation_landreman_paul_axi.f90)
target_link_libraries(plot_deviation_landreman_paul_axi.x myplot_lib
                                                          rabe_lib
                                                          helpers_lib)
add_test(PlotDeviationLandremanPaulAxi plot_deviation_landreman_paul_axi.x)
set_tests_properties(PlotDeviationLandremanPaulAxi PROPERTIES LABELS "external")

add_executable(plot_deviation_goodman_nautilus.x plot_deviation_goodman_nautilus.f90)
target_link_libraries(plot_deviation_goodman_nautilus.x myplot_lib
                                                          rabe_lib
                                                          helpers_lib)
add_test(PlotDeviationGoodmanNautilus plot_deviation_goodman_nautilus.x)
set_tests_properties(PlotDeviationGoodmanNautilus PROPERTIES LABELS "external")

add_executable(plot_deviation_drevlak_qh.x plot_deviation_drevlak_qh.f90)
target_link_libraries(plot_deviation_drevlak_qh.x myplot_lib
                                                          rabe_lib
                                                          helpers_lib)
add_test(PlotDeviationDrevlakQH plot_deviation_drevlak_qh.x)
set_tests_properties(PlotDeviationDrevlakQH PROPERTIES LABELS "external")

add_executable(plot_fieldlines_quasdex.x plot_fieldlines_quasdex.f90)
target_link_libraries(plot_fieldlines_quasdex.x myplot_lib
                                                rabe_lib
                                                helpers_lib)
add_test(PlotFieldlinesQuASDEX plot_fieldlines_quasdex.x)
set_tests_properties(PlotFieldlinesQuASDEX PROPERTIES LABELS "external")

add_executable(compare_distribution_function_neo2.x compare_distribution_function_neo2.f90)
target_link_libraries(compare_distribution_function_neo2.x myplot_lib
                                                          rabe_lib
                                                          helpers_lib)
add_test(CompareDistributionFunctionNeo2 compare_distribution_function_neo2.x)
set_tests_properties(CompareDistributionFunctionNeo2 PROPERTIES LABELS "external")

add_executable(compare_deviations.x compare_deviations.f90)
target_link_libraries(compare_deviations.x myplot_lib
                                           rabe_lib
                                           helpers_lib)
add_test(CompareDeviations compare_deviations.x)
set_tests_properties(CompareDeviations PROPERTIES LABELS "external")


add_executable(plot_helical_omnigenous.x plot_helical_omnigenous.f90)
target_link_libraries(plot_helical_omnigenous.x myplot_lib
                                                rabe_lib
                                                helpers_lib)
add_test(PlotHelicalOmnigenous plot_helical_omnigenous.x)
set_tests_properties(PlotHelicalOmnigenous PROPERTIES LABELS "external")

add_executable(plot_finite_collisionality_landreman_paul_axi.x plot_finite_collisionality_landreman_paul_axi.f90)
target_link_libraries(plot_finite_collisionality_landreman_paul_axi.x myplot_lib
                                                rabe_lib
                                                helpers_lib)
add_test(PlotFiniteCollisionality plot_finite_collisionality_landreman_paul_axi.x)
set_tests_properties(PlotFiniteCollisionality PROPERTIES LABELS "external")
