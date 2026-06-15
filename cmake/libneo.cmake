include(FetchContent)

# Override via -DLIBNEO_PATH=<dir> to use a local checkout instead of fetching.
# Override via -DLIBNEO_REF=<branch|tag|sha> to pin a specific libneo revision.
# Neither variable is read from the shell environment; set them only as CMake
# cache options to keep the build hermetic.

set(LIBNEO_PATH "" CACHE PATH "Local libneo source directory (leave empty to fetch)")
set(LIBNEO_REF  "" CACHE STRING "libneo git ref (branch, tag, or sha) to fetch")

if(LIBNEO_PATH AND EXISTS "${LIBNEO_PATH}")
    message(STATUS "Using libneo in ${LIBNEO_PATH}")
    add_subdirectory("${LIBNEO_PATH}"
                     "${CMAKE_CURRENT_BINARY_DIR}/libneo"
                     EXCLUDE_FROM_ALL)
else()
    if(LIBNEO_REF STREQUAL "")
        set(_libneo_ref "752aae9c8e31136121bbd77ce439f51e47a754be")
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
