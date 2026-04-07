add_library(precession_lib STATIC
  precession.f90
    bounce.f90
    orbit.f90
    params.f90
    field_instance.f90
  )

target_link_libraries(precession_lib PRIVATE utils_lib
field_lib
fieldline_lib
vmec_lib
LIBNEO::neo)
target_include_directories(precession_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
