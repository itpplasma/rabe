cmake_minimum_required(VERSION 3.24)

project(rabe)
enable_language(C Fortran)

set(CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
include(${CMAKE_DIR}/compiler.cmake)

add_subdirectory(src)
