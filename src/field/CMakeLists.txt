add_library(field_lib STATIC field_base.f90 fourier_field.f90)
target_link_libraries(field_lib PUBLIC utils_lib interpolate)
target_include_directories(field_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
