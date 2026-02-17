add_library(fieldline_lib STATIC
  fieldline.f90
  make_fieldline.f90
  fieldline_integrals.f90
  fieldline_integrands.f90
  surface_average.f90
  )

target_link_libraries(fieldline_lib PRIVATE utils_lib field_lib)
target_include_directories(fieldline_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
