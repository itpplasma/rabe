add_library(logger_lib STATIC
  logger.f90
  logger_config.f90
  )

target_include_directories(logger_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
