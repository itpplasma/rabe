program test_netcdf_output
    use constants, only: dp
    use netcdf_mod, only: netcdf_output_t, read_netcdf_values
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

    subroutine cleanup_test_file()
        logical :: file_exists
        integer :: unit, iostat

        inquire (file=test_file, exist=file_exists)
        if (file_exists) then
            open (newunit=unit, file=test_file, iostat=iostat)
            if (iostat == 0) then
                close (unit, status="delete")
            end if
        end if
    end subroutine cleanup_test_file

end program test_netcdf_output
