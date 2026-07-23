file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_against_simple.x test_against_simple.f90)
target_link_libraries(test_against_simple.x vmec_lib utils_lib)
add_test(TestAgainstSIMPLE test_against_simple.x)
set_tests_properties(TestAgainstSIMPLE PROPERTIES LABELS "quick")

add_executable(test_against_neo_qa.x test_against_neo_qa.f90)
target_link_libraries(test_against_neo_qa.x vmec_lib utils_lib)
add_test(TestAgainstNEO_QA test_against_neo_qa.x)
set_tests_properties(TestAgainstNEO_QA PROPERTIES LABELS "quick")

add_executable(test_against_neo_qh.x test_against_neo_qh.f90)
target_link_libraries(test_against_neo_qh.x vmec_lib utils_lib)
add_test(TestAgainstNEO_QH test_against_neo_qh.x)
set_tests_properties(TestAgainstNEO_QH PROPERTIES LABELS "quick")

add_executable(test_surface_area.x test_surface_area.f90)
target_link_libraries(test_surface_area.x vmec_lib utils_lib)
add_test(TestSurfaceArea test_surface_area.x)
set_tests_properties(TestSurfaceArea PROPERTIES LABELS "quick")

add_executable(test_stellarator_symmetry_qa.x test_stellarator_symmetry_qa.f90)
target_link_libraries(test_stellarator_symmetry_qa.x vmec_lib utils_lib)
add_test(TestStellaratorSymmetryQA test_stellarator_symmetry_qa.x)
set_tests_properties(TestStellaratorSymmetryQA PROPERTIES LABELS "quick")

add_executable(test_chartmap.x test_chartmap.f90)
target_link_libraries(test_chartmap.x vmec_lib utils_lib)
add_test(TestChartmap test_chartmap.x)
set_tests_properties(TestChartmap PROPERTIES LABELS "quick")

add_executable(test_booz_xform_signs.x test_booz_xform_signs.f90)
target_link_libraries(test_booz_xform_signs.x vmec_lib utils_lib LIBNEO::boozer)
add_test(TestBoozXformSigns test_booz_xform_signs.x)
set_tests_properties(TestBoozXformSigns PROPERTIES LABELS "quick")

# Diagnostic-only emitter consumed by plasma-sign-conventions.  It is not a
# standalone pass/fail test: the independent Cartesian verifier owns the gate.
add_executable(rabe_cartesian_sign_probe.x rabe_cartesian_sign_probe.f90)
target_link_libraries(rabe_cartesian_sign_probe.x vmec_lib field_lib fieldline_lib utils_lib)

# Diagnostic-only serializer emitter consumed by plasma-sign-conventions. The
# independent Python verifier decodes both files and owns the pass/fail gate.
add_executable(rabe_serializer_sign_probe.x rabe_serializer_sign_probe.f90)
target_link_libraries(
    rabe_serializer_sign_probe.x
    vmec_lib
    utils_lib
    netcdf_lib
    LIBNEO::boozer
)
