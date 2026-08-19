include(FetchContent)

set(FORTIO_REF "" CACHE STRING
    "fortio git ref (branch, tag, or SHA) to fetch")

if(NOT TARGET fortio)
    if(FORTIO_REF STREQUAL "")
        set(_fortio_ref "471256b7382afb8f641f37d833cc9172605da018")
    else()
        set(_fortio_ref "${FORTIO_REF}")
    endif()
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG ${_fortio_ref}
    )
    FetchContent_MakeAvailable(fortio)
endif()
