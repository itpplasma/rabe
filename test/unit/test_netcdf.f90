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
    character(len=*), parameter :: dim_name = "radius"
    integer, parameter :: n_dim = 3
    real(dp), parameter :: test_factor_a_array(n_dim) = [1.23_dp, 4.56_dp, 7.89_dp]
    real(dp), parameter :: test_factor_b_array(n_dim) = [98.7_dp, 6.54_dp, 32.1_dp]
    character(len=*), parameter :: title = "RABE Bootstrap Current Analysis Results"
    character(len=*), parameter :: long_name_a = "1/sqrt(nu_star) factor"
    character(len=*), parameter :: long_name_b = "1/nu_star factor"

    real(dp) :: read_factor_a, read_factor_b
    real(dp), dimension(n_dim) :: read_factor_a_array
    real(dp), dimension(n_dim) :: read_factor_b_array
    character(len=100) :: read_title, read_long_name_a, read_long_name_b

    logical :: file_exists
    integer :: unit, iostat

    call nc_out%create(test_file)
    call nc_out%add_global_attribute("title", title)
    call nc_out%add_real("off_factor_a")
    call nc_out%add_real_attr("off_factor_a", "long_name", long_name_a)
    call nc_out%add_real("off_factor_b")
    call nc_out%add_real_attr("off_factor_b", "long_name", long_name_b)
    call nc_out%write_real("off_factor_a", test_factor_a)
    call nc_out%write_real("off_factor_b", test_factor_b)
    call nc_out%close()

    inquire (file=test_file, exist=file_exists)
    if (.not. file_exists) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: file does not exist after creation"
        error stop
    end if

    call nc_in%open(test_file)
    call nc_in%read_real("off_factor_a", read_factor_a)
    call nc_in%read_real("off_factor_b", read_factor_b)
    call nc_in%read_global_attribute("title", read_title)
    call nc_in%read_real_attr("off_factor_a", "long_name", read_long_name_a)
    call nc_in%read_real_attr("off_factor_b", "long_name", read_long_name_b)
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

    if (trim(read_title) /= title) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: global title mismatch"
        print *, "found: ", trim(read_title)
        print *, "expected: RABE Bootstrap Current Analysis Results"
        error stop
    end if

    if (trim(read_long_name_a) /= long_name_a) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_a long_name mismatch"
        print *, "found: ", trim(read_long_name_a)
        print *, "expected: ", trim(long_name_a)
        error stop
    end if

    if (trim(read_long_name_b) /= long_name_b) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_b long_name mismatch"
        print *, "found: ", trim(read_long_name_b)
        print *, "expected: ", trim(long_name_b)
        error stop
    end if

    inquire (file=test_file, exist=file_exists)
    if (file_exists) then
        open (newunit=unit, file=test_file, iostat=iostat)
        if (iostat == 0) then
            close (unit, status="delete")
        end if
    else
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: file does not exist after reading"
        error stop
    end if

    call nc_out%create(test_file)
    call nc_out%add_global_attribute("title", title)
    call nc_out%def_dim(dim_name, n_dim)
    call nc_out%add_real_1d("off_factor_a", dim_name)
    call nc_out%add_real_attr("off_factor_a", "long_name", long_name_a)
    call nc_out%add_real_1d("off_factor_b", dim_name)
    call nc_out%add_real_attr("off_factor_b", "long_name", long_name_b)
    call nc_out%write_real_1d("off_factor_a", test_factor_a_array)
    call nc_out%write_real_1d("off_factor_b", test_factor_b_array)
    call nc_out%close()

    inquire (file=test_file, exist=file_exists)
    if (.not. file_exists) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: file does not exist after creation"
        error stop
    end if

    call nc_in%open(test_file)
    call nc_in%read_real_1d("off_factor_a", read_factor_a_array)
    call nc_in%read_real_1d("off_factor_b", read_factor_b_array)
    call nc_in%read_global_attribute("title", read_title)
    call nc_in%read_real_attr("off_factor_a", "long_name", read_long_name_a)
    call nc_in%read_real_attr("off_factor_b", "long_name", read_long_name_b)
    call nc_in%close()

    if (not_same(test_factor_a_array, read_factor_a_array, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_a_array mismatch"
        print *, "found: ", read_factor_a_array
        print *, "expected: ", test_factor_a_array
        error stop
    end if

    if (not_same(test_factor_b_array, read_factor_b_array, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_b mismatch"
        print *, "found: ", read_factor_b_array
        print *, "expected: ", test_factor_b_array
        error stop
    end if

    if (trim(read_title) /= title) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: global title mismatch"
        print *, "found: ", trim(read_title)
        print *, "expected: RABE Bootstrap Current Analysis Results"
        error stop
    end if

    if (trim(read_long_name_a) /= long_name_a) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_a long_name mismatch"
        print *, "found: ", trim(read_long_name_a)
        print *, "expected: ", trim(long_name_a)
        error stop
    end if

    if (trim(read_long_name_b) /= long_name_b) then
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: off_factor_b long_name mismatch"
        print *, "found: ", trim(read_long_name_b)
        print *, "expected: ", trim(long_name_b)
        error stop
    end if

    inquire (file=test_file, exist=file_exists)
    if (file_exists) then
        open (newunit=unit, file=test_file, iostat=iostat)
        if (iostat == 0) then
            close (unit, status="delete")
        end if
    else
        print *, "-------------------------------------------------------------"
        print *, "test_netcdf failed: file does not exist after reading"
        error stop
    end if

end program test_netcdf
