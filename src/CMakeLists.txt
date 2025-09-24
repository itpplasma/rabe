file(GLOB SRC_FILES "*.f90")

add_library(rabe_lib STATIC
  ${SRC_FILES}
  )

include(Util)

# Find NetCDF Fortran
find_package(PkgConfig REQUIRED)
pkg_check_modules(NetCDF_Fortran REQUIRED IMPORTED_TARGET netcdf-fortran)
target_include_directories(rabe_lib PRIVATE ${NetCDF_Fortran_INCLUDE_DIRS})
target_link_libraries(rabe_lib PRIVATE PkgConfig::NetCDF_Fortran)

find_or_fetch(quadpack)
target_compile_options(quadpack PRIVATE -w) # No warnings as library is not modified
target_link_libraries(rabe_lib PRIVATE quadpack)

add_subdirectory(neo)
target_link_libraries(rabe_lib PRIVATE neo_lib)

add_subdirectory(utils)
target_link_libraries(rabe_lib PRIVATE utils_lib)
