include(FetchContent)

if(DEFINED ENV{CODE} AND EXISTS "$ENV{CODE}/libneo")
    message(STATUS "Using libneo in $ENV{CODE}/libneo")
    add_subdirectory("$ENV{CODE}/libneo"
                     "${CMAKE_CURRENT_BINARY_DIR}/libneo"
                     EXCLUDE_FROM_ALL)
else()
    FetchContent_Declare(
        libneo
        GIT_REPOSITORY https://github.com/itpplasma/libneo.git
        GIT_TAG        main
        PATCH_COMMAND  ${CMAKE_COMMAND} -E copy
                       ${CMAKE_SOURCE_DIR}/cmake/libneo_CMakeLists.txt
                       <SOURCE_DIR>/CMakeLists.txt
    )
    FetchContent_MakeAvailable(libneo)
endif()
