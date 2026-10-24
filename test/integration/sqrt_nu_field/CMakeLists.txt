add_executable(test_sqrt_nu_field.x test_sqrt_nu_field.f90)
target_link_libraries(test_sqrt_nu_field.x sqrt_nu_geometry_lib helpers_lib rabe_lib)
add_test(NAME TestSqrtNuField COMMAND test_sqrt_nu_field.x)
set_tests_properties(TestSqrtNuField PROPERTIES LABELS "quick")
