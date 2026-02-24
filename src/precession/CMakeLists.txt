add_library(precession_lib STATIC
  precession.f90
  )

target_link_libraries(precession_lib PRIVATE utils_lib field_lib fieldline_lib)
target_include_directories(precession_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
