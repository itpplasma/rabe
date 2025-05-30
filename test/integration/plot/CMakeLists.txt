add_library(plot_lib STATIC
  plot_mod.f90
  )

add_subdirectory(pyplot-fortran)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/pyplot-fortran)
target_link_libraries(plot_lib PRIVATE pyplot_fortran)
