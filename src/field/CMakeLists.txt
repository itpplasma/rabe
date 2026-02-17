add_library(field_lib STATIC
  field_base.f90
  neo_field.f90
  )

target_link_libraries(field_lib PRIVATE utils_lib neo_lib)
target_include_directories(field_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
