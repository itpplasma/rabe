program test_netcdf
    use constants, only: dp
    use netcdf_mod, only: netcdf_t
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-12_dp

    type(netcdf_t) :: nc_out, nc_in
    character(len=*), parameter :: test_file = "test_output.nc"
    real(dp), parameter :: test_factor_a = 1.23456789_dp
    real(dp), parameter :: test_factor_b = 9.87654321_dp

    real(dp) :: read_factor_a, read_factor_b

    logical :: file_exists
    integer :: unit, iostat

    call nc_out%create(test_file)
    call nc_out%add_global_attribute("title", &
                                     "RABE Bootstrap Current Analysis Results")
    call nc_out%add_real("off_factor_a")
    call nc_out%add_real_attr("off_factor_a", "long_name", &
                              "1/sqrt(nu_star) factor")
    call nc_out%add_real("off_factor_b")
    call nc_out%add_real_attr("off_factor_b", "long_name", &
                              "1/nu_star factor")
    call nc_out%write_real("off_factor_a", test_factor_a)
    call nc_out%write_real("off_factor_b", test_factor_b)
    call nc_out%close()

    call nc_in%open(test_file)
    call nc_in%read_real("off_factor_a", read_factor_a)
    call nc_in%read_real("off_factor_b", read_factor_b)
    call nc_in%close()

    if (not_same(test_factor_a, read_factor_a, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_a mismatch"
        print *, "found: ", read_factor_a
        print *, "expected: ", test_factor_a
        error stop
    end if

    if (not_same(test_factor_b, read_factor_b, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_b mismatch"
        print *, "found: ", read_factor_b
        print *, "expected: ", test_factor_b
        error stop
    end if

    inquire (file=test_file, exist=file_exists)
    if (file_exists) then
        open (newunit=unit, file=test_file, iostat=iostat)
        if (iostat == 0) then
            close (unit, status="delete")
        end if
    end if

end program test_netcdf
