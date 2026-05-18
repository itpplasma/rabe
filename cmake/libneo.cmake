include(FetchContent)

if(DEFINED ENV{CODE} AND EXISTS "$ENV{CODE}/libneo")
    set(FETCHCONTENT_SOURCE_DIR_LIBNEO "$ENV{CODE}/libneo" CACHE PATH "")
    message(STATUS "Using libneo in $ENV{CODE}/libneo")
endif()

FetchContent_Declare(
    libneo
    GIT_REPOSITORY https://github.com/itpplasma/libneo.git
    GIT_TAG        main
    PATCH_COMMAND  ${CMAKE_COMMAND} -E copy
                   ${CMAKE_SOURCE_DIR}/cmake/libneo_CMakeLists.txt
                   <SOURCE_DIR>/CMakeLists.txt
)
FetchContent_MakeAvailable(libneo)
