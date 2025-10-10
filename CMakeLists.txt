cmake_minimum_required(VERSION 3.24)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(rabe)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(compiler)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/src/git_version.f90
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DOUTPUT_FILE=${CMAKE_BINARY_DIR}/src/git_version.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/git_version.cmake
    DEPENDS
        ${CMAKE_SOURCE_DIR}/src/git_version.f90.in
        ${CMAKE_SOURCE_DIR}/.git/HEAD
        ${CMAKE_SOURCE_DIR}/.git/index
    COMMENT "Generating git_version.f90 with current git hash"
    VERBATIM
)

add_custom_target(git_version_target
    DEPENDS ${CMAKE_BINARY_DIR}/src/git_version.f90
)

include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src/neo)
include_directories(${CMAKE_BINARY_DIR}/src/utils)
add_subdirectory(src)

add_executable(${PROJECT_NAME} app/main.f90)
target_link_libraries(${PROJECT_NAME} rabe_lib)

enable_testing()
add_subdirectory(test)
