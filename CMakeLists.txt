cmake_minimum_required(VERSION 3.24)

project(rabe)
enable_language(C Fortran)

set(CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
include(${CMAKE_DIR}/compiler.cmake)

set(SOURCES
    src/read_file.f90
    src/find_extrema.f90
)
add_library(rabe_lib ${SOURCES})

add_executable(${PROJECT_NAME} src/main.f90)
target_link_libraries(${PROJECT_NAME} rabe_lib)

enable_testing()
add_subdirectory(tests)
