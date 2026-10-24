add_executable(test_sqrt_nu_driver.x test_sqrt_nu_driver.f90)
target_link_libraries(test_sqrt_nu_driver.x rabe_lib)
add_test(NAME TestSqrtNuDriver COMMAND test_sqrt_nu_driver.x $<TARGET_FILE:sqrt_nu_rabe>)
set_tests_properties(TestSqrtNuDriver PROPERTIES LABELS "quick")
