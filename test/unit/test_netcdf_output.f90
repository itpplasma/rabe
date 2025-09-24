program test_netcdf_output
    use constants, only: dp
    use netcdf_output, only: netcdf_output_t, verify_netcdf_file

    implicit none

    type(netcdf_output_t) :: output
    character(len=*), parameter :: test_file = "test_output.nc"
    real(dp), parameter :: test_factor_a = 1.23456789_dp
    real(dp), parameter :: test_factor_b = 9.87654321_dp
    real(dp), parameter :: tolerance = 1.0e-12_dp

    logical :: verification_success

    print *, "Testing NetCDF output module..."

    call output%create(test_file)
    call output%write_results(test_factor_a, test_factor_b)
    call output%close()

    verification_success = verify_netcdf_file(test_file, test_factor_a, test_factor_b, tolerance)
    if (.not. verification_success) then
        error stop "Test failed"
    end if

    call cleanup_test_file()

    print *, "PASS: All NetCDF output tests passed"

contains

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