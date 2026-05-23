include(ExternalProject)
include(FindPkgConfig)

set(RABE_NETCDF_PROVIDER "auto" CACHE STRING "NetCDF provider: auto, system, or fetch")
set_property(CACHE RABE_NETCDF_PROVIDER PROPERTY STRINGS auto system fetch)

function(rabe_check_netcdf_fortran result_var)
    find_program(RABE_NF_CONFIG nf-config)
    if(NOT RABE_NF_CONFIG)
        set(${result_var} FALSE PARENT_SCOPE)
        return()
    endif()

    execute_process(
        COMMAND ${RABE_NF_CONFIG} --includedir
        OUTPUT_VARIABLE RABE_NETCDF_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    set(test_dir "${CMAKE_CURRENT_BINARY_DIR}/check-netcdf")
    set(test_src "${test_dir}/check_netcdf.f90")
    file(MAKE_DIRECTORY "${test_dir}")
    file(WRITE "${test_src}" "program check_netcdf\n  use netcdf\n  implicit none\nend program check_netcdf\n")

    execute_process(
        COMMAND ${CMAKE_Fortran_COMPILER} -I${RABE_NETCDF_INCLUDE_DIR} -fsyntax-only ${test_src}
        RESULT_VARIABLE check_result
        OUTPUT_QUIET
        ERROR_QUIET
        TIMEOUT 30
    )
    file(REMOVE_RECURSE "${test_dir}")

    if(check_result EQUAL 0)
        set(${result_var} TRUE PARENT_SCOPE)
        set(RABE_NETCDF_INCLUDE_DIR "${RABE_NETCDF_INCLUDE_DIR}" PARENT_SCOPE)
    else()
        set(${result_var} FALSE PARENT_SCOPE)
    endif()
endfunction()

function(rabe_add_fetched_netcdf)
    set(deps_prefix "${CMAKE_BINARY_DIR}/deps")
    set(deps_source_dir "${CMAKE_BINARY_DIR}/deps-src")
    set(deps_build_dir "${CMAKE_BINARY_DIR}/deps-build")
    set(netcdf_fortran_flags "${CMAKE_Fortran_FLAGS}")

    if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        string(APPEND netcdf_fortran_flags " -fPIE")
    endif()

    file(MAKE_DIRECTORY "${deps_prefix}/include")
    file(MAKE_DIRECTORY "${deps_prefix}/lib")

    set(common_cmake_args
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_Fortran_FLAGS=${netcdf_fortran_flags}
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_INSTALL_PREFIX=${deps_prefix}
        -DCMAKE_PREFIX_PATH=${deps_prefix}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DBUILD_SHARED_LIBS=OFF
    )

    ExternalProject_Add(rabe_hdf5_external
        URL https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
        URL_HASH SHA256=ec2e13c52e60f9a01491bb3158cb3778c985697131fc6a342262d32a26e58e44
        DOWNLOAD_DIR ${deps_source_dir}
        SOURCE_DIR ${deps_source_dir}/hdf5
        BINARY_DIR ${deps_build_dir}/hdf5
        INSTALL_DIR ${deps_prefix}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        CMAKE_ARGS
            ${common_cmake_args}
            -DHDF5_BUILD_FORTRAN=ON
            -DHDF5_BUILD_CPP_LIB=OFF
            -DHDF5_BUILD_JAVA=OFF
            -DHDF5_BUILD_EXAMPLES=OFF
            -DHDF5_BUILD_TOOLS=OFF
            -DHDF5_BUILD_HL_LIB=ON
            -DHDF5_ENABLE_PARALLEL=OFF
            -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
            -DBUILD_TESTING=OFF
            -DHDF5_BUILD_TESTING=OFF
        BUILD_BYPRODUCTS
            ${deps_prefix}/lib/libhdf5.a
            ${deps_prefix}/lib/libhdf5_hl.a
            ${deps_prefix}/lib/libhdf5_fortran.a
            ${deps_prefix}/lib/libhdf5_hl_fortran.a
            ${deps_prefix}/lib/libhdf5_f90cstub.a
            ${deps_prefix}/lib/libhdf5_hl_f90cstub.a
    )

    add_library(rabe_hdf5 STATIC IMPORTED GLOBAL)
    set_target_properties(rabe_hdf5 PROPERTIES
        IMPORTED_LOCATION ${deps_prefix}/lib/libhdf5.a
        INTERFACE_INCLUDE_DIRECTORIES ${deps_prefix}/include
        INTERFACE_LINK_LIBRARIES "z;dl;m"
    )
    add_dependencies(rabe_hdf5 rabe_hdf5_external)

    add_library(rabe_hdf5_hl STATIC IMPORTED GLOBAL)
    set_target_properties(rabe_hdf5_hl PROPERTIES
        IMPORTED_LOCATION ${deps_prefix}/lib/libhdf5_hl.a
        INTERFACE_INCLUDE_DIRECTORIES ${deps_prefix}/include
        INTERFACE_LINK_LIBRARIES "rabe_hdf5;z;dl;m"
    )
    add_dependencies(rabe_hdf5_hl rabe_hdf5_external)

    ExternalProject_Add(rabe_netcdf_c_external
        GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
        GIT_TAG v4.9.2
        GIT_SHALLOW TRUE
        DOWNLOAD_DIR ${deps_source_dir}
        SOURCE_DIR ${deps_source_dir}/netcdf-c
        BINARY_DIR ${deps_build_dir}/netcdf-c
        INSTALL_DIR ${deps_prefix}
        CMAKE_ARGS
            ${common_cmake_args}
            -DHDF5_ROOT=${deps_prefix}
            -DHDF5_DIR=${deps_prefix}
            -DENABLE_DAP=OFF
            -DENABLE_BYTERANGE=OFF
            -DENABLE_EXAMPLES=OFF
            -DENABLE_TESTS=OFF
            -DENABLE_EXTREME_NUMBERS=OFF
            -DENABLE_PARALLEL4=OFF
            -DENABLE_PNETCDF=OFF
            -DENABLE_CDF5=ON
            -DENABLE_NCZARR=OFF
            -DBUILD_UTILITIES=OFF
            -DENABLE_FILTER_TESTING=OFF
            -DENABLE_PLUGINS=OFF
            -DENABLE_HDF5=ON
            -DENABLE_NETCDF_4=ON
        DEPENDS rabe_hdf5_external
        BUILD_BYPRODUCTS ${deps_prefix}/lib/libnetcdf.a
    )

    add_library(netcdf::netcdf STATIC IMPORTED GLOBAL)
    set_target_properties(netcdf::netcdf PROPERTIES
        IMPORTED_LOCATION ${deps_prefix}/lib/libnetcdf.a
        INTERFACE_INCLUDE_DIRECTORIES ${deps_prefix}/include
        INTERFACE_LINK_LIBRARIES "rabe_hdf5_hl;rabe_hdf5;z;dl;m"
    )
    add_dependencies(netcdf::netcdf rabe_netcdf_c_external)

    ExternalProject_Add(rabe_netcdf_fortran_external
        GIT_REPOSITORY https://github.com/Unidata/netcdf-fortran.git
        GIT_TAG v4.6.1
        GIT_SHALLOW TRUE
        DOWNLOAD_DIR ${deps_source_dir}
        SOURCE_DIR ${deps_source_dir}/netcdf-fortran
        BINARY_DIR ${deps_build_dir}/netcdf-fortran
        INSTALL_DIR ${deps_prefix}
        CMAKE_ARGS
            ${common_cmake_args}
            -DNETCDF_C_LIBRARY=${deps_prefix}/lib/libnetcdf.a
            -DNETCDF_C_INCLUDE_DIR=${deps_prefix}/include
            -DENABLE_TESTS=OFF
            -DENABLE_EXAMPLES=OFF
        DEPENDS rabe_netcdf_c_external
        BUILD_BYPRODUCTS ${deps_prefix}/lib/libnetcdff.a
    )

    add_library(netcdf::netcdff STATIC IMPORTED GLOBAL)
    set_target_properties(netcdf::netcdff PROPERTIES
        IMPORTED_LOCATION ${deps_prefix}/lib/libnetcdff.a
        INTERFACE_INCLUDE_DIRECTORIES ${deps_prefix}/include
        INTERFACE_LINK_LIBRARIES netcdf::netcdf
    )
    add_dependencies(netcdf::netcdff rabe_netcdf_fortran_external)

    set(RABE_NETCDF_INCLUDE_DIR "${deps_prefix}/include" CACHE PATH "" FORCE)
endfunction()

function(rabe_configure_netcdf)
    if(RABE_NETCDF_PROVIDER STREQUAL "fetch")
        set(use_system FALSE)
    elseif(RABE_NETCDF_PROVIDER STREQUAL "system")
        set(use_system TRUE)
    elseif(RABE_NETCDF_PROVIDER STREQUAL "auto")
        rabe_check_netcdf_fortran(use_system)
    else()
        message(FATAL_ERROR "Invalid RABE_NETCDF_PROVIDER='${RABE_NETCDF_PROVIDER}'")
    endif()

    if(use_system)
        pkg_check_modules(NetCDF_Fortran REQUIRED IMPORTED_TARGET netcdf-fortran)
        if(NOT TARGET netcdf::netcdff)
            add_library(netcdf::netcdff INTERFACE IMPORTED GLOBAL)
            target_link_libraries(netcdf::netcdff INTERFACE PkgConfig::NetCDF_Fortran)
        endif()
        if(NOT RABE_NETCDF_INCLUDE_DIR)
            find_program(RABE_NF_CONFIG nf-config)
            execute_process(
                COMMAND ${RABE_NF_CONFIG} --includedir
                OUTPUT_VARIABLE RABE_NETCDF_INCLUDE_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        endif()
        message(STATUS "Using system NetCDF-Fortran")
    else()
        message(STATUS "Building NetCDF-Fortran for ${CMAKE_Fortran_COMPILER_ID}")
        rabe_add_fetched_netcdf()
    endif()
endfunction()
