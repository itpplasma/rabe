cmake_minimum_required(VERSION 3.24)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(rabe)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(compiler)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(SuiteSparse)
if(SUITESPARSE_FOUND)
  # In Debian bookworm being standard in 2025, SuiteSparse v5.12 packages don't have
  # their own CMake files for find_package. As soon as this changes,
  # find_package(SuiteSparse) should be replaced by find_package(UMFPACK) and
  # the legacy switch removed.
  set(LEGACY_SUITESPARSE ON)
  message(STATUS "SuiteSparse found.")
else()
  # Use more recent SuiteSparse v7.6.1 from GitHub.
  message(STATUS "SuiteSparse not found. Building.")
  include(suitesparse)
endif()

include_directories(${CMAKE_BINARY_DIR}/src)
add_subdirectory(src)

add_executable(${PROJECT_NAME} app/main.f90)
target_link_libraries(${PROJECT_NAME} rabe_lib)

enable_testing()
add_subdirectory(test)
