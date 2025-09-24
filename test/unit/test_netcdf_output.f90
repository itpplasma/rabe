program test_netcdf_output
    use constants, only: dp
    use netcdf_output, only: netcdf_output_t
    use netcdf
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-12_dp

    type(netcdf_output_t) :: output
    character(len=*), parameter :: test_file = "test_output.nc"
    real(dp), parameter :: test_factor_a = 1.23456789_dp
    real(dp), parameter :: test_factor_b = 9.87654321_dp

    real(dp) :: read_factor_a, read_factor_b

    call output%create(test_file)
    call output%write_results(test_factor_a, test_factor_b)
    call output%close()

    call read_netcdf_values(test_file, read_factor_a, read_factor_b)

    if (not_same(test_factor_a, read_factor_a, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf_output failed: off_factor_a mismatch"
        print *, "found: ", read_factor_a
        print *, "expected: ", test_factor_a
        error stop
    end if

    if (not_same(test_factor_b, read_factor_b, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf_output failed: off_factor_b mismatch"
        print *, "found: ", read_factor_b
        print *, "expected: ", test_factor_b
        error stop
    end if

    call cleanup_test_file()

contains

    subroutine read_netcdf_values(filename, factor_a, factor_b)
        character(len=*), intent(in) :: filename
        real(dp), intent(out) :: factor_a, factor_b

        integer :: ncid, var_id_a, var_id_b, status

        status = nf90_open(filename, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            print *, "-------------------------------------------------------------"
            print *, "test_netcdf_output failed: Cannot open NetCDF file"
            error stop
        end if

        status = nf90_inq_varid(ncid, "off_factor_a", var_id_a)
        if (status /= NF90_NOERR) then
            print *, "-------------------------------------------------------------"
            print *, "test_netcdf_output failed: Variable off_factor_a not found"
            status = nf90_close(ncid)
            error stop
        end if

        status = nf90_get_var(ncid, var_id_a, factor_a)
        if (status /= NF90_NOERR) then
            print *, "-------------------------------------------------------------"
            print *, "test_netcdf_output failed: Cannot read off_factor_a"
            status = nf90_close(ncid)
            error stop
        end if

        status = nf90_inq_varid(ncid, "off_factor_b", var_id_b)
        if (status /= NF90_NOERR) then
            print *, "-------------------------------------------------------------"
            print *, "test_netcdf_output failed: Variable off_factor_b not found"
            status = nf90_close(ncid)
            error stop
        end if

        status = nf90_get_var(ncid, var_id_b, factor_b)
        if (status /= NF90_NOERR) then
            print *, "-------------------------------------------------------------"
            print *, "test_netcdf_output failed: Cannot read off_factor_b"
            status = nf90_close(ncid)
            error stop
        end if

        status = nf90_close(ncid)
        if (status /= NF90_NOERR) then
            print *, "WARNING: Error closing NetCDF file"
        end if
    end subroutine read_netcdf_values

    subroutine cleanup_test_file()
        logical :: file_exists
        integer :: unit, iostat

        inquire(file=test_file, exist=file_exists)
        if (file_exists) then
            open(newunit=unit, file=test_file, iostat=iostat)
            if (iostat == 0) then
                close(unit, status="delete")
            end if
        end if
    end subroutine cleanup_test_file

end program test_netcdf_output