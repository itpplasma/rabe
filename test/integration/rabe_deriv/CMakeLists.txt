# rabe_deriv.x on a field, then on its own output as spectrum file, and on a
# wout with the libneo VMEC -> Boozer chain (wout_derivative).
set(RD_DIR ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../vmec/input/wout_LandremanPaul2021_QH_reactorScale_ns101_reference.nc
     DESTINATION ${RD_DIR})
configure_file(rabe_from_field.in ${RD_DIR}/from_field/rabe.in COPYONLY)
configure_file(rabe_from_spectrum.in ${RD_DIR}/from_spectrum/rabe.in COPYONLY)
configure_file(rabe_from_wout.in ${RD_DIR}/from_wout/rabe.in COPYONLY)

add_executable(compare_rabe_deriv.x compare_rabe_deriv.f90)
target_link_libraries(compare_rabe_deriv.x netcdf_lib utils_lib)

add_test(NAME RabeDerivFromField COMMAND ${CMAKE_BINARY_DIR}/rabe_deriv.x
         WORKING_DIRECTORY ${RD_DIR}/from_field)
set_tests_properties(RabeDerivFromField PROPERTIES LABELS "quick"
                     FIXTURES_SETUP RabeDerivField)
add_test(NAME RabeDerivFromSpectrum COMMAND ${CMAKE_BINARY_DIR}/rabe_deriv.x
         WORKING_DIRECTORY ${RD_DIR}/from_spectrum)
set_tests_properties(RabeDerivFromSpectrum PROPERTIES LABELS "quick"
                     FIXTURES_REQUIRED RabeDerivField
                     FIXTURES_SETUP RabeDerivSpectrum)
add_test(NAME RabeDerivRoundTrip COMMAND compare_rabe_deriv.x
         WORKING_DIRECTORY ${RD_DIR})
set_tests_properties(RabeDerivRoundTrip PROPERTIES LABELS "quick"
                     FIXTURES_REQUIRED "RabeDerivField;RabeDerivSpectrum")

add_test(NAME RabeDerivFromWout COMMAND ${CMAKE_BINARY_DIR}/rabe_deriv.x
         WORKING_DIRECTORY ${RD_DIR}/from_wout)
set_tests_properties(RabeDerivFromWout PROPERTIES LABELS "quick")
