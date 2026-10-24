cmake_minimum_required(VERSION 3.24)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(rabe)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(compiler)

# required for fetching dependencies
find_package(Git REQUIRED)

include(fortio)

include(git_version_target)

add_subdirectory(src) # for rabe_lib

add_executable(${PROJECT_NAME} app/main.f90)
target_link_libraries(${PROJECT_NAME} rabe_lib)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "rabe.x")

enable_testing()
add_subdirectory(test)

option(BUILD_PYTHON_BINDINGS "Build Python bindings via f90wrap" OFF)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()
