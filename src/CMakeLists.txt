set(SOURCES
    read_file.f90
    find_extrema.f90
)
add_library(rabe_lib ${SOURCES})

add_executable(${PROJECT_NAME} main.f90)
target_link_libraries(${PROJECT_NAME} rabe_lib)