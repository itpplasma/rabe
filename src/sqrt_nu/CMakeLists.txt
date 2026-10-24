add_library(sqrt_nu_lib STATIC
    sqrt_nu.f90
    ../generated/sqrt_nu_kernel.f90
)
target_include_directories(sqrt_nu_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
find_package(OpenMP COMPONENTS Fortran QUIET)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(sqrt_nu_lib PRIVATE OpenMP::OpenMP_Fortran)
endif()
