add_library(netcdf_lib STATIC
    netcdf.f90
  )

target_link_libraries(netcdf_lib PUBLIC utils_lib)
target_link_libraries(netcdf_lib PUBLIC netcdf::netcdff)
target_include_directories(netcdf_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
