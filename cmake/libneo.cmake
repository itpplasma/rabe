include(FetchContent)

# -DLIBNEO_PATH=<dir>: build against a local libneo checkout instead of fetching.
set(LIBNEO_PATH "" CACHE PATH "Local libneo source directory (leave empty to fetch)")

# -DLIBNEO_REF=<branch|tag|sha>: fetch this libneo revision instead of the pin.
set(LIBNEO_REF  "" CACHE STRING "libneo git ref (branch, tag, or sha) to fetch")

# rabe links only libneo's light targets (boozer, vmec_support, interpolate);
# skip its numerics-heavy core so BLAS/LAPACK and fortnum are not required here.
set(LIBNEO_BUILD_NUMERICS OFF CACHE BOOL "rabe needs only light libneo targets" FORCE)

if(LIBNEO_PATH AND EXISTS "${LIBNEO_PATH}")
    message(STATUS "Using libneo in ${LIBNEO_PATH}")
    add_subdirectory("${LIBNEO_PATH}"
                     "${CMAKE_CURRENT_BINARY_DIR}/libneo"
                     EXCLUDE_FROM_ALL)
else()
    if(LIBNEO_REF STREQUAL "")
        set(_libneo_ref "94fd53bfc3e45420d48b830978674efbe37d8dcc")
    else()
        set(_libneo_ref "${LIBNEO_REF}")
    endif()
    FetchContent_Declare(
        libneo
        GIT_REPOSITORY https://github.com/itpplasma/libneo.git
        GIT_TAG        ${_libneo_ref}
        EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(libneo)
endif()
