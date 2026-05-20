include(FetchContent)

if(DEFINED ENV{LIBNEO} AND EXISTS "$ENV{LIBNEO}")
    message(STATUS "Using libneo in $ENV{LIBNEO}")
    add_subdirectory("$ENV{LIBNEO}"
                     "${CMAKE_CURRENT_BINARY_DIR}/libneo"
                     EXCLUDE_FROM_ALL)
else()
    FetchContent_Declare(
        libneo
        GIT_REPOSITORY https://github.com/itpplasma/libneo.git
        GIT_TAG        752aae9c8e31136121bbd77ce439f51e47a754be
        PATCH_COMMAND  ${CMAKE_COMMAND} -E copy
                       ${CMAKE_SOURCE_DIR}/cmake/libneo_CMakeLists.txt
                       <SOURCE_DIR>/CMakeLists.txt
    )
    FetchContent_MakeAvailable(libneo)
endif()
