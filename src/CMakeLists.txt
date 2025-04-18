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

add_library(rabe_lib STATIC
  ${SRC_FILES}
  ${SUITESPARSE_SRC_FILES}
  )

target_link_libraries(rabe_lib PUBLIC
BLAS::BLAS
LAPACK::LAPACK
ZLIB::ZLIB
Threads::Threads
)

if(LEGACY_SUITESPARSE)
  target_link_libraries(rabe_lib PUBLIC ${SUITESPARSE_LIBRARIES})
else()
  target_link_libraries(rabe_lib PUBLIC SuiteSparse::umfpack)
endif()
