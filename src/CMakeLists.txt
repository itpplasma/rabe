file(GLOB SRC_FILES "*.f90")

add_library(rabe_lib STATIC
  ${SRC_FILES}
  ${CMAKE_BINARY_DIR}/src/git_version.f90
  )

add_dependencies(rabe_lib git_version_target)
target_include_directories(rabe_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Find NetCDF Fortran
find_package(PkgConfig REQUIRED)
pkg_check_modules(NetCDF_Fortran REQUIRED IMPORTED_TARGET netcdf-fortran)
target_include_directories(rabe_lib PRIVATE ${NetCDF_Fortran_INCLUDE_DIRS})
target_link_libraries(rabe_lib PRIVATE PkgConfig::NetCDF_Fortran)

add_subdirectory(neo)
add_subdirectory(utils)
add_subdirectory(field)
add_subdirectory(fieldline)
add_subdirectory(shaing_callen)

target_link_libraries(rabe_lib PUBLIC
  utils_lib neo_lib field_lib fieldline_lib shaing_callen_lib)
