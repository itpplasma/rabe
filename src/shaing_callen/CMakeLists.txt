add_library(shaing_callen_lib STATIC
  shaing_callen.f90
  shaing_callen_integration.f90
  shaing_callen_wrappers.f90
  )

target_link_libraries(shaing_callen_lib PRIVATE
  utils_lib field_lib fieldline_lib)
target_include_directories(shaing_callen_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
