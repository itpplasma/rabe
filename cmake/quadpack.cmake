include(FetchContent)

FetchContent_Declare(
    quadpack
    GIT_REPOSITORY https://github.com/jacobwilliams/quadpack.git
    GIT_TAG        702abfd5f0acbdb51439695334347a4b3c0dc87a
    PATCH_COMMAND  ${CMAKE_COMMAND} -E copy
                   ${CMAKE_CURRENT_LIST_DIR}/quadpack_CMakeLists.txt
                   <SOURCE_DIR>/CMakeLists.txt
)
FetchContent_MakeAvailable(quadpack)
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-w" RABE_FORTRAN_SUPPORTS_LOWERCASE_W)
if(RABE_FORTRAN_SUPPORTS_LOWERCASE_W)
    target_compile_options(quadpack PRIVATE -w)
endif()
