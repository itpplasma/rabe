include(FetchContent)

# -DFORTIO_REF=<branch|tag|SHA>: fetch this Fortio revision instead of the pin.
set(FORTIO_REF "829fcde02024047f05e7a2389a01709617525e13" CACHE STRING
    "fortio git ref (branch, tag, or SHA) to fetch")

if(NOT TARGET fortio)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG ${FORTIO_REF}
    )
    FetchContent_MakeAvailable(fortio)
endif()
