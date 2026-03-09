add_library(vmec_lib STATIC
            boozer_converter.F90
            magfie.f90
            magfie_can_boozer.f90
            field_vmec.f90
            )

find_or_fetch(libneo)
target_link_libraries(vmec_lib PRIVATE utils_lib LIBNEO::neo)
target_include_directories(vmec_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
