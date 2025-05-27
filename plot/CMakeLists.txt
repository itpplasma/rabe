cmake_minimum_required (VERSION 3.20)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(rabe-plot)

add_subdirectory(pyplot-fortran)
add_executable(plot_maxima plot_maxima.f90)
target_link_libraries(plot_maxima pyplot_fortran rabe_lib)

add_custom_target(plot
DEPENDS plot_maxima
)
