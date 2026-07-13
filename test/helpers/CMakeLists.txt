add_library(helpers_lib STATIC
  mock_field.f90
  mock_perturbed_field.f90
  anti_sigma_field.f90
  horner_fourier_field.f90
  plot_quantities.f90
  readers.f90
  test_shaing_callen_mod.f90
  bc_file.f90
)
# myplot_lib PUBLIC: flang (unlike gfortran) requires myplot_module.mod (and
# transitively pyplot_module.mod) to be resolvable wherever plot_quantities.mod
# is used, since plot_quantities.f90 uses myplot_module in its public interface.
target_link_libraries(helpers_lib PRIVATE offset_lib
                                          utils_lib
                                          field_lib
                                          fourier_field_lib
                                          vmec_lib
                                          fieldline_lib
                                          shaing_callen_lib
                                  PUBLIC  myplot_lib)
