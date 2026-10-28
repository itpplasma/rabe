add_library(derivative_lib STATIC
  spectral_surface.f90
  spectrum_derivative.f90
  spectrum_io.f90
  wout_chain.f90
  )

# wout_chain exposes libneo's boozer_vmec_ad types, so LIBNEO::boozer is public.
target_link_libraries(derivative_lib PRIVATE utils_lib offset_lib netcdf_lib
  LIBNEO::vmec_support
  PUBLIC field_lib fieldline_lib LIBNEO::boozer)
target_include_directories(derivative_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
