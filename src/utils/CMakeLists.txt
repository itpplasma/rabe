add_library(utils_lib STATIC
  utils.f90
  constants.f90
  integrate.f90
  fourier.f90
  find_extrema.f90
  diophantine.f90
  trace.f90
  )

include(quadpack)
target_link_libraries(utils_lib PRIVATE quadpack)
target_include_directories(utils_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
