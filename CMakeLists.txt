cmake_minimum_required(VERSION 3.24)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(rabe)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(compiler)

# OFF (default): normal build - netcdf_lib/vmec_lib and app/main.f90 (rabe.x) as usual.
# ON: NetCDF-free AD build - only the Fourier pipeline (utils/field/fieldline/
# offset/shaing_callen/error_handling + libneo's interpolate module, see
# cmake/libneo.cmake), no VMEC/NetCDF/rabe.x/tests, plus python/ad -> _rabe_ad.
# Intended for a flang toolchain (no prebuilt netcdf-fortran, AD path never
# needs VMEC). Not meant to be set directly by users of the normal build - see
# RABE_WITH_AD_BACKEND, whose sub-build sets this automatically.
option(RABE_AD_BUILD "NetCDF-free AD build: only the Fourier pipeline + _rabe_ad" OFF)

# required for fetching dependencies
find_package(Git REQUIRED)

if(EXISTS ${PROJECT_SOURCE_DIR}/.git/HEAD)
    set(GIT_HEAD_DEPS
        ${PROJECT_SOURCE_DIR}/.git/HEAD
        ${PROJECT_SOURCE_DIR}/.git/index
    )
endif()

add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/src/git_version.f90
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DOUTPUT_FILE=${PROJECT_BINARY_DIR}/src/git_version.f90
        -P ${PROJECT_SOURCE_DIR}/cmake/git_version.cmake
    DEPENDS
        ${PROJECT_SOURCE_DIR}/src/git_version.f90.in
        ${GIT_HEAD_DEPS}
    COMMENT "Generating git_version.f90 with current git hash"
    VERBATIM
)

add_custom_target(git_version_target
    DEPENDS ${PROJECT_BINARY_DIR}/src/git_version.f90
)

add_subdirectory(src) # for rabe_lib

if(NOT RABE_AD_BUILD)
    add_executable(${PROJECT_NAME} app/main.f90)
    target_link_libraries(${PROJECT_NAME} rabe_lib)
    set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "rabe.x")

    enable_testing()
    add_subdirectory(test)
endif()

option(BUILD_PYTHON_BINDINGS "Build Python bindings via f90wrap" OFF)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(RABE_AD_BUILD)
    add_subdirectory(python/ad)
endif()
