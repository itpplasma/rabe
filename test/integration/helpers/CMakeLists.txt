add_library(helpers_lib STATIC
  mock_field.f90
  mock_perturbed_field.f90
  anti_sigma_field.f90
  plot_quantities.f90
)
target_link_libraries(helpers_lib PRIVATE rabe_lib
                                          myplot_lib)
