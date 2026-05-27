add_library(rabe_lib STATIC
  rabe.f90
  ${CMAKE_CURRENT_BINARY_DIR}/git_version.f90
)

add_dependencies(rabe_lib git_version_target)
target_include_directories(rabe_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_library(read_file_lib STATIC read_file.f90)
target_include_directories(read_file_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(read_file_lib PRIVATE utils_lib fieldline_lib)

add_subdirectory(error_handling)
add_subdirectory(utils)
add_subdirectory(netcdf)
add_subdirectory(vmec)
add_subdirectory(field)
add_subdirectory(fieldline)
add_subdirectory(offset)
add_subdirectory(shaing_callen)

target_link_libraries(rabe_lib PRIVATE
  error_handling_lib
  utils_lib
  netcdf_lib
  field_lib
  vmec_lib
  fieldline_lib
  shaing_callen_lib
  offset_lib)
