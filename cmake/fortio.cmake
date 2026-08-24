include(FetchContent)

if(NOT TARGET fortio)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG 829fcde02024047f05e7a2389a01709617525e13
    )
    FetchContent_MakeAvailable(fortio)
endif()
