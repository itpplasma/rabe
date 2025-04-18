file(GLOB SRC_FILES "*.f90")

add_library(rabe_lib STATIC
  ${SRC_FILES}
  )

add_subdirectory(neo)
target_link_libraries(rabe_lib PRIVATE neo_lib)
