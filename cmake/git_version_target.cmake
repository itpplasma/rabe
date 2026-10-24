# Git may store HEAD in a worktree gitfile and the active ref in loose or
# packed refs. Check the actual revision on each build instead of enumerating
# those mutable paths. configure_file writes only when the hash has changed.
add_custom_target(git_version_target
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DOUTPUT_FILE=${PROJECT_BINARY_DIR}/src/git_version.f90
        -P ${CMAKE_CURRENT_LIST_DIR}/git_version.cmake
    BYPRODUCTS ${PROJECT_BINARY_DIR}/src/git_version.f90
    COMMENT "Checking git_version.f90 against current git hash"
    VERBATIM
)
