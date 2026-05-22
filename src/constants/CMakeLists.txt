add_library(constants_lib STATIC constants.f90)
target_include_directories(constants_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
