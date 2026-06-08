include(FetchContent)

if(DEFINED ENV{LIBNEO_PATH} AND EXISTS "$ENV{LIBNEO_PATH}")
    message(STATUS "Using libneo in $ENV{LIBNEO_PATH}")
    add_subdirectory("$ENV{LIBNEO_PATH}"
                     "${CMAKE_CURRENT_BINARY_DIR}/libneo"
                     EXCLUDE_FROM_ALL)
else()
    # LIBNEO_REF env overrides the pinned ref so a release can test a candidate.
    if(DEFINED ENV{LIBNEO_REF} AND NOT "$ENV{LIBNEO_REF}" STREQUAL "")
        set(_libneo_ref "$ENV{LIBNEO_REF}")
    else()
        set(_libneo_ref "752aae9c8e31136121bbd77ce439f51e47a754be")
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
