add_library(utils_lib STATIC
  utils.f90
  constants.f90
  integrate.f90
  fourier.f90
  find_extrema.f90
  diophantine.f90
  )

include(quadpack)
# PUBLIC: quadpack_double is a private implementation detail of integrate.f90
# (qag is renamed on import and never appears in integrate's public interface),
# but flang (unlike gfortran) requires quadpack_double.mod to be transitively
# resolvable wherever `integrate` is used, not just within utils_lib itself
# (e.g. fieldline_integrals.f90).
target_link_libraries(utils_lib PUBLIC quadpack)
target_include_directories(utils_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
