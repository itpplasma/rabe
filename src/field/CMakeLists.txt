add_library(field_lib STATIC field_base.f90)
target_link_libraries(field_lib PRIVATE constants_lib logger_lib)
target_include_directories(field_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
