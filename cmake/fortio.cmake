include(FetchContent)

if(NOT TARGET fortio)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG 4af4f714eecc1430879c95dd144f6ead229158ce
    )
    FetchContent_MakeAvailable(fortio)
endif()
