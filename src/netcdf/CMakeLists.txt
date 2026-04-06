add_library(netcdf_lib STATIC
    netcdf.f90
  )

# Find NetCDF Fortran
find_package(PkgConfig REQUIRED)
pkg_check_modules(NetCDF_Fortran REQUIRED IMPORTED_TARGET netcdf-fortran)
target_include_directories(netcdf_lib PUBLIC ${NetCDF_Fortran_INCLUDE_DIRS})

target_link_libraries(netcdf_lib PRIVATE utils_lib PkgConfig::NetCDF_Fortran)
target_include_directories(netcdf_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
