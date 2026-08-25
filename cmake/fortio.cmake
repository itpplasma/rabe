include(FetchContent)

# -DFORTIO_REF=<branch|tag|SHA>: fetch this Fortio revision instead of the pin.
set(FORTIO_REF "4af4f714eecc1430879c95dd144f6ead229158ce" CACHE STRING
    "fortio git ref (branch, tag, or SHA) to fetch")

if(NOT TARGET fortio)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG ${FORTIO_REF}
    )
    FetchContent_MakeAvailable(fortio)
endif()
