add_library(vmec_lib STATIC
            boozer_field.f90
            )

include(libneo)
target_link_libraries(vmec_lib PRIVATE utils_lib field_lib
                      LIBNEO::boozer LIBNEO::vmec_field
                      LIBNEO::interpolate LIBNEO::vmec_support)
target_include_directories(vmec_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
