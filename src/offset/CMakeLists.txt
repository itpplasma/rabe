add_library(offset_lib STATIC
  coefficients.f90
  deviation.f90
  fit_functions.f90
)

target_include_directories(offset_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(offset_lib PRIVATE
  utils_lib
  error_handling_lib
  field_lib
  fieldline_lib
  logger_lib)
