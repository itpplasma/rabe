module netcdf_output
    use constants, only: dp
    use netcdf

    implicit none
    private

    public :: netcdf_output_t, verify_netcdf_file

    type :: netcdf_output_t
        integer :: ncid = -1
        logical :: is_open = .false.
        integer :: var_id_a = -1
        integer :: var_id_b = -1
    contains
        procedure :: create => netcdf_output_create
        procedure :: write_results => netcdf_output_write_results
        procedure :: close => netcdf_output_close
        final :: netcdf_output_final
    end type netcdf_output_t

contains

    subroutine netcdf_output_create(this, filename)
        class(netcdf_output_t), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: status

        if (this%is_open) then
            call this%close()
        end if

        status = nf90_create(filename, NF90_CLASSIC_MODEL, this%ncid)
        call check_netcdf_status(status, "creating file: " // filename)

        status = nf90_def_var(this%ncid, "off_factor_a", NF90_DOUBLE, &
                              varid=this%var_id_a)
        call check_netcdf_status(status, "defining variable off_factor_a")

        status = nf90_put_att(this%ncid, this%var_id_a, "long_name", &
                              "1/sqrt(nu_star) factor")
        call check_netcdf_status(status, "setting attribute for off_factor_a")

        status = nf90_def_var(this%ncid, "off_factor_b", NF90_DOUBLE, &
                              varid=this%var_id_b)
        call check_netcdf_status(status, "defining variable off_factor_b")

        status = nf90_put_att(this%ncid, this%var_id_b, "long_name", &
                              "1/nu_star factor")
        call check_netcdf_status(status, "setting attribute for off_factor_b")

        status = nf90_put_att(this%ncid, NF90_GLOBAL, "title", &
                              "RABE Bootstrap Current Analysis Results")
        call check_netcdf_status(status, "setting global title")

        status = nf90_enddef(this%ncid)
        call check_netcdf_status(status, "ending definition mode")

        this%is_open = .true.
    end subroutine netcdf_output_create

    subroutine netcdf_output_write_results(this, off_factor_a, off_factor_b)
        class(netcdf_output_t), intent(inout) :: this
        real(dp), intent(in) :: off_factor_a, off_factor_b
        integer :: status

        if (.not. this%is_open) then
            error stop "NetCDF file not open for writing"
        end if

        status = nf90_put_var(this%ncid, this%var_id_a, off_factor_a)
        call check_netcdf_status(status, "writing off_factor_a")

        status = nf90_put_var(this%ncid, this%var_id_b, off_factor_b)
        call check_netcdf_status(status, "writing off_factor_b")

    end subroutine netcdf_output_write_results

    subroutine netcdf_output_close(this)
        class(netcdf_output_t), intent(inout) :: this
        integer :: status

        if (this%is_open) then
            status = nf90_close(this%ncid)
            call check_netcdf_status(status, "closing NetCDF file")
            this%is_open = .false.
            this%ncid = -1
        end if
    end subroutine netcdf_output_close

    subroutine netcdf_output_final(this)
        type(netcdf_output_t), intent(inout) :: this

        if (this%is_open) then
            call this%close()
        end if
    end subroutine netcdf_output_final

    subroutine check_netcdf_status(status, operation)
        integer, intent(in) :: status
        character(len=*), intent(in) :: operation

        if (status /= NF90_NOERR) then
            print *, "NetCDF error during ", operation, ":"
            print *, trim(nf90_strerror(status))
            error stop "NetCDF operation failed"
        end if
    end subroutine check_netcdf_status

    function verify_netcdf_file(filename, expected_a, expected_b, tolerance) result(success)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: expected_a, expected_b, tolerance
        logical :: success

        integer :: ncid, var_id_a, var_id_b, status
        real(dp) :: read_a, read_b
        logical :: file_exists

        success = .false.

        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: NetCDF file does not exist: ", filename
            return
        end if

        status = nf90_open(filename, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            print *, "FAIL: Cannot open NetCDF file: ", filename
            return
        end if

        status = nf90_inq_varid(ncid, "off_factor_a", var_id_a)
        if (status /= NF90_NOERR) then
            print *, "FAIL: Variable off_factor_a not found"
            status = nf90_close(ncid)
            return
        end if

        status = nf90_get_var(ncid, var_id_a, read_a)
        if (status /= NF90_NOERR) then
            print *, "FAIL: Cannot read off_factor_a"
            status = nf90_close(ncid)
            return
        end if

        status = nf90_inq_varid(ncid, "off_factor_b", var_id_b)
        if (status /= NF90_NOERR) then
            print *, "FAIL: Variable off_factor_b not found"
            status = nf90_close(ncid)
            return
        end if

        status = nf90_get_var(ncid, var_id_b, read_b)
        if (status /= NF90_NOERR) then
            print *, "FAIL: Cannot read off_factor_b"
            status = nf90_close(ncid)
            return
        end if

        status = nf90_close(ncid)
        if (status /= NF90_NOERR) then
            print *, "WARNING: Error closing NetCDF file"
        end if

        if (abs(read_a - expected_a) > tolerance) then
            print *, "FAIL: off_factor_a mismatch"
            print *, "Expected:", expected_a, "Got:", read_a
            return
        end if

        if (abs(read_b - expected_b) > tolerance) then
            print *, "FAIL: off_factor_b mismatch"
            print *, "Expected:", expected_b, "Got:", read_b
            return
        end if

        success = .true.
    end function verify_netcdf_file

end module netcdf_output