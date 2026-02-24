file(GLOB SRC_FILES "*.f90")

add_library(rabe_lib STATIC
  ${SRC_FILES}
  ${CMAKE_BINARY_DIR}/src/git_version.f90
  )

add_dependencies(rabe_lib git_version_target)
target_include_directories(rabe_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(neo)
add_subdirectory(utils)
add_subdirectory(netcdf)
add_subdirectory(vmec)
add_subdirectory(field)
add_subdirectory(fieldline)
add_subdirectory(shaing_callen)
add_subdirectory(precession)

target_link_libraries(rabe_lib PUBLIC
  utils_lib
  neo_lib
  netcdf_lib
  field_lib
  vmec_lib
  fieldline_lib
  shaing_callen_lib
  precession_lib)
