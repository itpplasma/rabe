include(FetchContent)

# -DFORTIO_REF=<branch|tag|SHA>: fetch this Fortio revision instead of the pin.
set(FORTIO_REF "471256b7382afb8f641f37d833cc9172605da018" CACHE STRING
    "fortio git ref (branch, tag, or SHA) to fetch")

if(NOT TARGET fortio)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG ${FORTIO_REF}
    )
    FetchContent_MakeAvailable(fortio)
endif()
