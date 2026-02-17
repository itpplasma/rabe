find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(SuiteSparse)
if(SUITESPARSE_FOUND)
  # In Debian bookworm being standard in 2025, SuiteSparse v5.12 packages don't have
  # their own CMake files for find_package. As soon as this changes,
  # find_package(SuiteSparse) should be replaced by find_package(UMFPACK) and
  # the legacy switch removed.
  set(LEGACY_SUITESPARSE ON)
  message(STATUS "SuiteSparse found.")
else()
  # Use more recent SuiteSparse v7.6.1 from GitHub.
  message(STATUS "SuiteSparse not found. Building.")
  include(suitesparse)
endif()

file(GLOB SRC_FILES "*.f90")

set(SUITESPARSE_SRC_FILES
  umf4_f77wrapper_ver_4_5.c
  umf4_f77zwrapper_ver_4_5.c
)
if(LEGACY_SUITESPARSE)
  set(SUITESPARSE_FLAGS "-I${SUITESPARSE_INCLUDE_DIRS}")
  message(${SUITESPARSE_INCLUDE_DIRS})
else()
  set(SUITESPARSE_FLAGS "")
endif()
set_source_files_properties(umf4_f77wrapper_ver_4_5.c  PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DDLONG")
set_source_files_properties(umf4_f77zwrapper_ver_4_5.c PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DZLONG")

add_library(neo_lib STATIC
  ${SRC_FILES}
  ${SUITESPARSE_SRC_FILES}
  )

target_link_libraries(neo_lib PRIVATE
BLAS::BLAS
LAPACK::LAPACK
ZLIB::ZLIB
Threads::Threads
)

if(LEGACY_SUITESPARSE)
  target_link_libraries(neo_lib PRIVATE ${SUITESPARSE_LIBRARIES})
else()
  target_link_libraries(neo_lib PRIVATE SuiteSparse::umfpack)
endif()

target_compile_options(neo_lib PRIVATE -w) # No warnings as legacy code is not modified
target_include_directories(neo_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
