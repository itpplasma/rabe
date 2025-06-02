add_subdirectory(plot)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/plot)

enable_testing()
add_subdirectory(unit)
add_subdirectory(integration)
