# Builds _rabe_ad: a plain shared library exposing the flat bind(C) AD driver
# (ad_driver.f90) over the NetCDF-free core (utils/field/fieldline/offset).
# No f90wrap involved - loaded from Python via ctypes (see python/ad/rabe_ad).
# Intended to be built with RABE_AD_BUILD=ON and a flang compiler
# (see cmake/ad_backend.cmake), but also builds fine with gfortran for local
# testing of the driver itself.

add_library(_rabe_ad SHARED ad_driver.f90)

target_link_libraries(_rabe_ad PRIVATE
    utils_lib
    error_handling_lib
    field_lib
    fourier_field_lib
    interpolate
    fieldline_lib
    offset_lib
)

set_target_properties(_rabe_ad PROPERTIES
    PREFIX ""
    OUTPUT_NAME "_rabe_ad"
)

install(TARGETS _rabe_ad DESTINATION rabe_ad)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/rabe_ad/__init__.py DESTINATION rabe_ad)
