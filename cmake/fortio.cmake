include(FetchContent)

if(NOT TARGET fortio)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG 11afd0bd1af0c99ea4e9a1c0df683dbdcc299b69
    )
    FetchContent_MakeAvailable(fortio)
endif()
