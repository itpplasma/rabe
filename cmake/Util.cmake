include(FetchContent)

function(find_or_fetch DEPENDENCY)
    cmake_parse_arguments(OPTIONAL "" "SUBDIR" "" ${ARGN})

    if(DEFINED ENV{CODE} AND EXISTS $ENV{CODE}/${DEPENDENCY})
        set(SOURCE_DIR $ENV{CODE}/${DEPENDENCY}/${OPTIONAL_SUBDIR})
        message(STATUS "Using ${DEPENDENCY}/${OPTIONAL_SUBDIR} in $ENV{CODE}/${DEPENDENCY}")
    else()
        fetch(${DEPENDENCY} SOURCE_DIR SUBDIR ${OPTIONAL_SUBDIR})
    endif()

    add_subdirectory(${SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/${DEPENDENCY}/${OPTIONAL_SUBDIR}
        EXCLUDE_FROM_ALL
    )
endfunction()

function(fetch DEPENDENCY SOURCE_DIR)
    cmake_parse_arguments(OPTIONAL "" "SUBDIR" "" ${ARGN})

    set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)
    get_branch_or_main(${REPO_URL} REMOTE_BRANCH)
    message(STATUS "Fetch ${DEPENDENCY} branch ${REMOTE_BRANCH} from ${REPO_URL}")


    string(TOLOWER ${DEPENDENCY} DEPENDENCY_LOWER)

    FetchContent_GetProperties(${DEPENDENCY})
    if(NOT ${DEPENDENCY_LOWER}_POPULATED)
        FetchContent_Declare(
            ${DEPENDENCY}
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            GIT_REPOSITORY ${REPO_URL}
            GIT_TAG ${REMOTE_BRANCH}
        )
        FetchContent_Populate(${DEPENDENCY})
    endif()
    message(STATUS "Fetched ${DEPENDENCY} to ${${DEPENDENCY_LOWER}_SOURCE_DIR}")

    set(${SOURCE_DIR} ${${DEPENDENCY_LOWER}_SOURCE_DIR}/${OPTIONAL_SUBDIR} PARENT_SCOPE)
endfunction()

function(get_branch_or_main REPO_URL REMOTE_BRANCH)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND git ls-remote --heads ${REPO_URL} ${BRANCH}
        OUTPUT_VARIABLE BRANCH_EXISTS
        ERROR_VARIABLE GIT_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message(STATUS "Branch ${BRANCH} exists: ${BRANCH_EXISTS}")

    if(BRANCH_EXISTS)
        set(${REMOTE_BRANCH} ${BRANCH} PARENT_SCOPE)
    else()
        set(${REMOTE_BRANCH} "main" PARENT_SCOPE)
    endif()
endfunction()
