add_library(helpers_lib STATIC
  mock_field.f90
  mock_perturbed_field.f90
  anti_sigma_field.f90
  plot_quantities.f90
  readers.f90
  distribution_function.f90
  test_cumint_mod.f90
  test_shaing_callen_mod.f90
)
target_link_libraries(helpers_lib PRIVATE rabe_lib
                                          myplot_lib)
