program test_netcdf_output
    use constants, only: dp
    use netcdf_output, only: netcdf_output_t
#ifdef HAVE_NETCDF
    use netcdf
#endif

    implicit none

    type(netcdf_output_t) :: output
    character(len=*), parameter :: test_file = "test_output.nc"
    real(dp), parameter :: test_factor_a = 1.23456789_dp
    real(dp), parameter :: test_factor_b = 9.87654321_dp
    real(dp), parameter :: tolerance = 1.0e-12_dp

#ifdef HAVE_NETCDF
    real(dp) :: read_factor_a, read_factor_b
    integer :: ncid, var_id, status
    logical :: file_exists
#endif

    print *, "Testing NetCDF output module..."

#ifdef HAVE_NETCDF

    call output%create(test_file)
    call output%write_results(test_factor_a, test_factor_b)
    call output%close()

    inquire(file=test_file, exist=file_exists)
    if (.not. file_exists) then
        print *, "FAIL: NetCDF file was not created"
        error stop "Test failed"
    end if

    status = nf90_open(test_file, NF90_NOWRITE, ncid)
    if (status /= NF90_NOERR) then
        print *, "FAIL: Cannot open test file for reading"
        print *, trim(nf90_strerror(status))
        error stop "Test failed"
    end if

    status = nf90_inq_varid(ncid, "off_factor_a", var_id)
    if (status /= NF90_NOERR) then
        print *, "FAIL: Variable off_factor_a not found"
        error stop "Test failed"
    end if

    status = nf90_get_var(ncid, var_id, read_factor_a)
    if (status /= NF90_NOERR) then
        print *, "FAIL: Cannot read off_factor_a"
        error stop "Test failed"
    end if

    status = nf90_inq_varid(ncid, "off_factor_b", var_id)
    if (status /= NF90_NOERR) then
        print *, "FAIL: Variable off_factor_b not found"
        error stop "Test failed"
    end if

    status = nf90_get_var(ncid, var_id, read_factor_b)
    if (status /= NF90_NOERR) then
        print *, "FAIL: Cannot read off_factor_b"
        error stop "Test failed"
    end if

    status = nf90_close(ncid)
    if (status /= NF90_NOERR) then
        print *, "WARNING: Cannot close test file"
    end if

    if (abs(read_factor_a - test_factor_a) > tolerance) then
        print *, "FAIL: off_factor_a mismatch"
        print *, "Expected:", test_factor_a
        print *, "Got:", read_factor_a
        error stop "Test failed"
    end if

    if (abs(read_factor_b - test_factor_b) > tolerance) then
        print *, "FAIL: off_factor_b mismatch"
        print *, "Expected:", test_factor_b
        print *, "Got:", read_factor_b
        error stop "Test failed"
    end if

    call cleanup_test_file()

    print *, "PASS: All NetCDF output tests passed"
#else
    print *, "SKIP: NetCDF support not available - testing stub functionality"
    call output%create(test_file)
    call output%write_results(test_factor_a, test_factor_b)
    call output%close()
    print *, "PASS: NetCDF stub functionality works"
#endif

contains

    subroutine cleanup_test_file()
#ifdef HAVE_NETCDF
        integer :: unit, iostat
        open(newunit=unit, file=test_file, iostat=iostat)
        if (iostat == 0) then
            close(unit, status="delete")
        end if
#endif
    end subroutine cleanup_test_file

end program test_netcdf_output