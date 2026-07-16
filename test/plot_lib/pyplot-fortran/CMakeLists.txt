add_library(pyplot_fortran STATIC
pyplot_module.F90
)
target_include_directories(pyplot_fortran PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
