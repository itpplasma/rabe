add_library(helpers_lib STATIC
  mock_field.f90
)
target_link_libraries(helpers_lib PRIVATE rabe_lib)
