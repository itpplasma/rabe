add_subdirectory(plot_lib)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/plot_lib)

add_subdirectory(helpers)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/helpers)

enable_testing()
add_subdirectory(unit)
add_subdirectory(integration)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_subdirectory(plot)
else()
    message(STATUS "Python3 not found - plot tests disabled")
endif()
add_subdirectory(golden)
