add_library(numerics_lib STATIC
  integrate.f90
  fourier.f90
  find_extrema.f90
  )

include(Util)
find_or_fetch(quadpack)
target_compile_options(quadpack PRIVATE -w)
target_link_libraries(numerics_lib PRIVATE utils_lib quadpack)
target_include_directories(numerics_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
