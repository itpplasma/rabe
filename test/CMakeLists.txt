add_subdirectory(plot)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/plot)

add_subdirectory(helpers)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/helpers)

enable_testing()
add_subdirectory(unit)
add_subdirectory(integration)
add_subdirectory(external)
