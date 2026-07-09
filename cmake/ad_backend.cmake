include(ExternalProject)

# Builds the flang/Enzyme AD backend (_rabe_ad) as an out-of-tree sub-build
# with its own Fortran compiler. This is required, not optional: one CMake
# configure has exactly one CMAKE_Fortran_COMPILER, so the flang-compiled
# NetCDF-free extension cannot be produced in the same configure as the
# gfortran-compiled f90wrap extension (_rabe). The two meet only at the
# Python `import` level once both are installed side by side.
#
# The sub-build re-configures this same source tree with RABE_AD_BUILD=ON
# (skips netcdf_lib/vmec_lib/rabe.x/tests, builds libneo's interpolate module
# in a NetCDF-free way - see cmake/libneo.cmake's RABE_AD_BUILD branch - and
# builds python/ad -> _rabe_ad.so, see python/ad/CMakeLists.txt). Its install
# step stages results locally; the outer install() below copies that staging
# tree into the wheel/package install destination (${SKBUILD_PROJECT_NAME}/ad),
# so this only takes effect once `cmake --install` runs, same as any other
# target in this project. Included from python/CMakeLists.txt, so
# SKBUILD_PROJECT_NAME is already defined by the time this runs.

if(NOT DEFINED RABE_AD_FORTRAN_COMPILER)
    set(RABE_AD_FORTRAN_COMPILER "flang-new" CACHE STRING
        "Fortran compiler (LLVM flang) used to build the rabe.ad AD backend")
endif()

set(RABE_AD_BACKEND_STAGE_DIR "${CMAKE_BINARY_DIR}/ad_backend_install")

ExternalProject_Add(rabe_ad_backend
    SOURCE_DIR      "${PROJECT_SOURCE_DIR}"
    BINARY_DIR      "${CMAKE_BINARY_DIR}/ad_backend_build"
    CMAKE_ARGS
        -DCMAKE_Fortran_COMPILER=${RABE_AD_FORTRAN_COMPILER}
        -DCMAKE_BUILD_TYPE=Release
        -DRABE_AD_BUILD=ON
        -DBUILD_PYTHON_BINDINGS=OFF
        -DRABE_WITH_AD_BACKEND=OFF
        -DCMAKE_INSTALL_PREFIX=${RABE_AD_BACKEND_STAGE_DIR}
    BUILD_COMMAND   ${CMAKE_COMMAND} --build <BINARY_DIR> --target _rabe_ad
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>
    BUILD_ALWAYS    TRUE
)

# Make the sub-build part of the normal build graph (so `cmake --build .`,
# as invoked by `pip install .`, actually produces it before `cmake --install`
# runs on the outer project).
add_custom_target(rabe_ad_backend_all ALL DEPENDS rabe_ad_backend)

install(DIRECTORY "${RABE_AD_BACKEND_STAGE_DIR}/rabe_ad/"
        DESTINATION "${SKBUILD_PROJECT_NAME}/ad"
        FILES_MATCHING PATTERN "*.so" PATTERN "*.py")
