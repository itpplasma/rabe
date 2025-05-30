file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_against_neo2_field.x test_against_neo2_field.f90)
target_link_libraries(test_against_neo2_field.x rabe_lib)
add_test(TestAgainstNeo2Field test_against_neo2_field.x)
set_tests_properties(TestAgainstNeo2Field PROPERTIES LABELS "slow")

add_executable(test_fieldline.x test_fieldline.f90)
target_link_libraries(test_fieldline.x rabe_lib)
add_test(TestFieldline test_fieldline.x)

add_executable(test_make_flock_of_fieldlines.x test_make_flock_of_fieldlines.f90)
target_link_libraries(test_make_flock_of_fieldlines.x rabe_lib)
add_test(TestMakeFlockOfFieldlines test_make_flock_of_fieldlines.x)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/helpers)
add_subdirectory(helpers)

add_executable(test_integration_along_fieldline.x test_integration_along_fieldline.f90)
target_link_libraries(test_integration_along_fieldline.x rabe_lib helpers_lib)
add_test(TestIntegrationFieldline test_integration_along_fieldline.x)

add_executable(test_find_global_maximum.x test_find_global_maximum.f90)
target_link_libraries(test_find_global_maximum.x rabe_lib helpers_lib)
add_test(TestGlobalMaximum test_find_global_maximum.x)

add_executable(test_anti_sigma_analytic.x test_anti_sigma_analytic.f90)
target_link_libraries(test_anti_sigma_analytic.x rabe_lib helpers_lib)
add_test(TestAntiSigmaAnalytic test_anti_sigma_analytic.x)

add_executable(test_pert_anti_sigma_analytic.x test_pert_anti_sigma_analytic.f90)
target_link_libraries(test_pert_anti_sigma_analytic.x rabe_lib helpers_lib)
add_test(TestPertAntiSigmaAnalytic test_pert_anti_sigma_analytic.x)

add_subdirectory(plot)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/plot)

add_executable(test_plot.x test_plot.f90)
target_link_libraries(test_plot.x myplot_lib utils_lib)
add_test(TestPlot test_plot.x)
set_tests_properties(TestPlot PROPERTIES LABELS "per_hand")

add_executable(test_plot_anti_sigma.x test_plot_anti_sigma.f90)
target_link_libraries(test_plot_anti_sigma.x myplot_lib
                                             rabe_lib
                                             helpers_lib)
add_test(TestPlotAntiSigma test_plot_anti_sigma.x)
set_tests_properties(TestPlotAntiSigma PROPERTIES LABELS "plot")
