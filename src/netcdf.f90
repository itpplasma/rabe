module netcdf_mod
    use constants, only: dp
    use netcdf

    implicit none
    private

    public :: netcdf_output_t, read_netcdf_values

    integer, parameter :: MAX_VARS = 100

    type :: var_info_t
        character(len=100) :: name
        integer :: var_id
    end type var_info_t

    type :: netcdf_output_t
        integer :: ncid = -1
        logical :: is_open = .false.
        logical :: in_define_mode = .false.
        integer :: n_vars = 0
        type(var_info_t) :: vars(MAX_VARS)
    contains
        procedure :: create => netcdf_output_create
        procedure :: add_global_attribute => netcdf_output_add_global_attr
        procedure :: add_real => netcdf_output_add_real
        procedure :: add_real_attr => netcdf_output_add_real_attr
        procedure :: end_define => netcdf_output_end_define
        procedure :: write_real => netcdf_output_write_real
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
        call check_netcdf_status(status, "creating file: "//filename)

        this%is_open = .true.
        this%in_define_mode = .true.
        this%n_vars = 0
    end subroutine netcdf_output_create

    subroutine netcdf_output_add_global_attr(this, attr_name, attr_value)
        class(netcdf_output_t), intent(inout) :: this
        character(len=*), intent(in) :: attr_name
        character(len=*), intent(in) :: attr_value
        integer :: status

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        status = nf90_put_att(this%ncid, NF90_GLOBAL, attr_name, attr_value)
        call check_netcdf_status(status, "setting global attribute: " &
            //attr_name)
    end subroutine netcdf_output_add_global_attr

    subroutine netcdf_output_add_real(this, var_name)
        class(netcdf_output_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        integer :: status, var_id

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        if (this%n_vars >= MAX_VARS) then
            error stop "Maximum number of variables exceeded"
        end if

        status = nf90_def_var(this%ncid, var_name, NF90_DOUBLE, varid=var_id)
        call check_netcdf_status(status, "defining variable: "//var_name)

        this%n_vars = this%n_vars + 1
        this%vars(this%n_vars)%name = var_name
        this%vars(this%n_vars)%var_id = var_id
    end subroutine netcdf_output_add_real

    subroutine netcdf_output_add_real_attr(this, var_name, attr_name, &
        attr_value)
        class(netcdf_output_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        character(len=*), intent(in) :: attr_name
        character(len=*), intent(in) :: attr_value
        integer :: status, var_id, i

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            error stop "Not in define mode"
        end if

        var_id = -1
        do i = 1, this%n_vars
            if (trim(this%vars(i)%name) == trim(var_name)) then
                var_id = this%vars(i)%var_id
                exit
            end if
        end do

        if (var_id == -1) then
            error stop "Variable not found: "//var_name
        end if

        status = nf90_put_att(this%ncid, var_id, attr_name, attr_value)
        call check_netcdf_status(status, "setting attribute "//attr_name &
            //" for variable: "//var_name)
    end subroutine netcdf_output_add_real_attr

    subroutine netcdf_output_end_define(this)
        class(netcdf_output_t), intent(inout) :: this
        integer :: status

        if (.not. this%is_open) then
            error stop "NetCDF file not open"
        end if

        if (.not. this%in_define_mode) then
            return
        end if

        status = nf90_enddef(this%ncid)
        call check_netcdf_status(status, "ending definition mode")

        this%in_define_mode = .false.
    end subroutine netcdf_output_end_define

    subroutine netcdf_output_write_real(this, var_name, value)
        class(netcdf_output_t), intent(inout) :: this
        character(len=*), intent(in) :: var_name
        real(dp), intent(in) :: value
        integer :: status, var_id, i

        if (.not. this%is_open) then
            error stop "NetCDF file not open for writing"
        end if

        if (this%in_define_mode) then
            call this%end_define()
        end if

        var_id = -1
        do i = 1, this%n_vars
            if (trim(this%vars(i)%name) == trim(var_name)) then
                var_id = this%vars(i)%var_id
                exit
            end if
        end do

        if (var_id == -1) then
            error stop "Variable not found: "//var_name
        end if

        status = nf90_put_var(this%ncid, var_id, value)
        call check_netcdf_status(status, "writing variable: "//var_name)
    end subroutine netcdf_output_write_real

    subroutine netcdf_output_close(this)
        class(netcdf_output_t), intent(inout) :: this
        integer :: status

        if (this%is_open) then
            if (this%in_define_mode) then
                call this%end_define()
            end if

            status = nf90_close(this%ncid)
            call check_netcdf_status(status, "closing NetCDF file")
            this%is_open = .false.
            this%in_define_mode = .false.
            this%ncid = -1
            this%n_vars = 0
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

    subroutine read_netcdf_values(filename, factor_a, factor_b)
        character(len=*), intent(in) :: filename
        real(dp), intent(out) :: factor_a, factor_b

        integer :: ncid, var_id_a, var_id_b, status

        status = nf90_open(filename, NF90_NOWRITE, ncid)
        call check_netcdf_status(status, "opening file: "//filename)

        status = nf90_inq_varid(ncid, "off_factor_a", var_id_a)
        call check_netcdf_status(status, "finding variable off_factor_a")

        status = nf90_get_var(ncid, var_id_a, factor_a)
        call check_netcdf_status(status, "reading off_factor_a")

        status = nf90_inq_varid(ncid, "off_factor_b", var_id_b)
        call check_netcdf_status(status, "finding variable off_factor_b")

        status = nf90_get_var(ncid, var_id_b, factor_b)
        call check_netcdf_status(status, "reading off_factor_b")

        status = nf90_close(ncid)
        call check_netcdf_status(status, "closing NetCDF file")
    end subroutine read_netcdf_values

end module netcdf_mod
