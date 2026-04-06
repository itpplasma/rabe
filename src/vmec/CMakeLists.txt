add_library(vmec_lib STATIC
            boozer_converter.F90
            boozer_field.f90
            )

find_or_fetch(libneo)
target_link_libraries(vmec_lib PRIVATE utils_lib field_lib LIBNEO::neo)
target_include_directories(vmec_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
