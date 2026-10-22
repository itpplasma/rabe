add_library(derivative_lib STATIC
  spectral_surface.f90
  spectrum_derivative.f90
  spectrum_io.f90
  )

target_link_libraries(derivative_lib PRIVATE utils_lib offset_lib netcdf_lib
  PUBLIC field_lib fieldline_lib)
target_include_directories(derivative_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
