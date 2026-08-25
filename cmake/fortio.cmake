include(FetchContent)

# -DFORTIO_REF=<branch|tag|SHA>: fetch this Fortio revision instead of the pin.
set(FORTIO_REF "11afd0bd1af0c99ea4e9a1c0df683dbdcc299b69" CACHE STRING
    "fortio git ref (branch, tag, or SHA) to fetch")

if(NOT TARGET fortio)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG ${FORTIO_REF}
    )
    FetchContent_MakeAvailable(fortio)
endif()
