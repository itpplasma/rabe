add_executable(test_sqrt_nu.x test_sqrt_nu.f90)
target_link_libraries(test_sqrt_nu.x PRIVATE sqrt_nu_lib)
add_test(NAME TestSqrtNuExact COMMAND test_sqrt_nu.x)
set_tests_properties(TestSqrtNuExact PROPERTIES LABELS "quick" ENVIRONMENT "OMP_NUM_THREADS=2")
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(NAME TestSqrtNuModelOracles
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/log_bvp_reference.py
                $<TARGET_FILE:test_sqrt_nu.x>)
set_tests_properties(TestSqrtNuModelOracles PROPERTIES
        LABELS "quick" ENVIRONMENT "OMP_NUM_THREADS=2")
