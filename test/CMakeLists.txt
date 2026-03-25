add_subdirectory(plot_lib)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/plot_lib)

add_subdirectory(helpers)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/helpers)

enable_testing()
add_subdirectory(unit)
add_subdirectory(integration)
add_subdirectory(external)
add_subdirectory(plot)
add_subdirectory(golden)
