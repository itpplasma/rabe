add_executable(test_sd22_f1_rabe_sign.x test_sd22_f1_rabe_sign.f90)
target_link_libraries(test_sd22_f1_rabe_sign.x helpers_lib utils_lib fieldline_lib
                      shaing_callen_lib rabe_lib)
add_test(NAME TestSD22F1RabeSign COMMAND test_sd22_f1_rabe_sign.x)
set_tests_properties(TestSD22F1RabeSign PROPERTIES LABELS "quick")

# NEO side of F1 needs a NEO executable (not part of rabe); run on request with
#   cmake -DSD22_F1_NEO_EXECUTABLE=/path/neo.x -DSD22_F1_NEO_SOURCE=/path/NEO ...
if(SD22_F1_NEO_EXECUTABLE AND SD22_F1_NEO_SOURCE)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_test(NAME TestSD22F1NeoSign
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/neo_f1_sign.py
                     ${SD22_F1_NEO_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/neo_f1
                     ${SD22_F1_NEO_SOURCE})
    set_tests_properties(TestSD22F1NeoSign PROPERTIES LABELS "neo")
endif()
