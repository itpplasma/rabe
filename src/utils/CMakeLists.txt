add_library(utils_lib STATIC
  utils.f90
  constants.f90
  )

target_include_directories(utils_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
