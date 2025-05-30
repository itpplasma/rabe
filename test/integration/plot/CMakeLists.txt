add_library(myplot_lib STATIC
  myplot_module.f90
  )

add_subdirectory(pyplot-fortran)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/pyplot-fortran)
target_link_libraries(myplot_lib PRIVATE pyplot_fortran)
target_link_libraries(myplot_lib PRIVATE utils_lib)
