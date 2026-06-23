include(FetchContent)

# -DLIBNEO_PATH=<dir>: build against a local libneo checkout instead of fetching.
set(LIBNEO_PATH "" CACHE PATH "Local libneo source directory (leave empty to fetch)")

# -DLIBNEO_REF=<branch|tag|sha>: fetch this libneo revision instead of the pin.
set(LIBNEO_REF  "" CACHE STRING "libneo git ref (branch, tag, or sha) to fetch")

if(LIBNEO_PATH AND EXISTS "${LIBNEO_PATH}")
    message(STATUS "Using libneo in ${LIBNEO_PATH}")
    add_subdirectory("${LIBNEO_PATH}"
                     "${CMAKE_CURRENT_BINARY_DIR}/libneo"
                     EXCLUDE_FROM_ALL)
else()
    if(LIBNEO_REF STREQUAL "")
        set(_libneo_ref "34b1b33f0705da84fef2d42482e4f96eb2a364d6")
    else()
        set(_libneo_ref "${LIBNEO_REF}")
    endif()
    FetchContent_Declare(
        libneo
        GIT_REPOSITORY https://github.com/itpplasma/libneo.git
        GIT_TAG        ${_libneo_ref}
        PATCH_COMMAND  ${CMAKE_COMMAND} -E copy
                       ${CMAKE_CURRENT_LIST_DIR}/libneo_CMakeLists.txt
                       <SOURCE_DIR>/CMakeLists.txt
    )
    FetchContent_MakeAvailable(libneo)
endif()
