module netcdf_output
    use constants, only: dp
#ifdef HAVE_NETCDF
    use netcdf
#endif

    implicit none
    private

    public :: netcdf_output_t

    type :: netcdf_output_t
        integer :: ncid = -1
        logical :: is_open = .false.
#ifdef HAVE_NETCDF
        integer :: var_id_a = -1
        integer :: var_id_b = -1
#endif
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
#ifdef HAVE_NETCDF
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
#else
        print *, "Warning: NetCDF support not available, cannot create ", filename
        print *, "         Consider installing NetCDF Fortran library"
#endif
    end subroutine netcdf_output_create

    subroutine netcdf_output_write_results(this, off_factor_a, off_factor_b)
        class(netcdf_output_t), intent(inout) :: this
        real(dp), intent(in) :: off_factor_a, off_factor_b

#ifdef HAVE_NETCDF
        integer :: status

        if (.not. this%is_open) then
            error stop "NetCDF file not open for writing"
        end if

        status = nf90_put_var(this%ncid, this%var_id_a, off_factor_a)
        call check_netcdf_status(status, "writing off_factor_a")

        status = nf90_put_var(this%ncid, this%var_id_b, off_factor_b)
        call check_netcdf_status(status, "writing off_factor_b")
#else
        print *, "Warning: NetCDF support not available, cannot write results"
        print *, "         Values: off_factor_a=", off_factor_a, " off_factor_b=", off_factor_b
#endif

    end subroutine netcdf_output_write_results

    subroutine netcdf_output_close(this)
        class(netcdf_output_t), intent(inout) :: this
#ifdef HAVE_NETCDF
        integer :: status

        if (this%is_open) then
            status = nf90_close(this%ncid)
            call check_netcdf_status(status, "closing NetCDF file")
            this%is_open = .false.
            this%ncid = -1
        end if
#endif
    end subroutine netcdf_output_close

    subroutine netcdf_output_final(this)
        type(netcdf_output_t), intent(inout) :: this

        if (this%is_open) then
            call this%close()
        end if
    end subroutine netcdf_output_final

#ifdef HAVE_NETCDF
    subroutine check_netcdf_status(status, operation)
        integer, intent(in) :: status
        character(len=*), intent(in) :: operation

        if (status /= NF90_NOERR) then
            print *, "NetCDF error during ", operation, ":"
            print *, trim(nf90_strerror(status))
            error stop "NetCDF operation failed"
        end if
    end subroutine check_netcdf_status
#endif

end module netcdf_output