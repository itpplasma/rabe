add_library(field_lib STATIC field_base.f90)
target_link_libraries(field_lib PRIVATE utils_lib)
target_include_directories(field_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_library(fourier_field_lib STATIC fourier_field.f90)
target_link_libraries(fourier_field_lib PRIVATE field_lib interpolate utils_lib)
target_include_directories(fourier_field_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
