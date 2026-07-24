add_library(vmec_lib STATIC
            boozer_field.f90
            )

include(libneo)
# PUBLIC: boozer_field_t (and callers like app/main.f90 that declare it
# directly) need libneo's boozer/vmec_support modules transitively resolvable
# under flang.
target_link_libraries(vmec_lib PUBLIC utils_lib field_lib
                      LIBNEO::boozer
                      LIBNEO::vmec_support)
target_include_directories(vmec_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
