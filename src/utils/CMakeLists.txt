add_library(utils_lib STATIC
  utils.f90
  constants.f90
  integrate.f90
  fourier.f90
  find_extrema.f90
  diophantine.f90
  )

include(quadpack)
# quadpack_double is a private implementation detail of integrate.f90 (qag is
# renamed on import and never appears in integrate's public interface), so
# PRIVATE is the correct visibility here - unlike src/field's interpolate link,
# nothing in utils_lib's own public API exposes a quadpack type. flang, unlike
# gfortran, still requires quadpack_double.mod to be transitively resolvable
# wherever `integrate` is used, so the AD build (RABE_AD_BUILD) needs it PUBLIC
# as a flang-tooling workaround; the normal (gfortran) build keeps the correct
# PRIVATE encapsulation.
if(RABE_AD_BUILD)
    target_link_libraries(utils_lib PUBLIC quadpack)
else()
    target_link_libraries(utils_lib PRIVATE quadpack)
endif()
target_include_directories(utils_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
