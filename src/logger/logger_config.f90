module logger_config

    implicit none

    character(len=256), public, protected :: log_file
    character(len=256), public, protected :: log_level

    namelist /log_config/ &
        log_file, &
        log_level

contains

    subroutine read_namelist(filename)
        character(len=*), intent(in) :: filename
        integer :: ios, unit
        logical :: file_exists

        !> Defaults
        log_file = "rabe.log"
        log_level = "INFO"

        inquire (file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "Error in logger_config: file not found: ", filename
            error stop
        end if

        open (newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print *, "Error in logger_config: cannot open file: ", filename
            error stop
        end if

        read (unit, nml=log_config, iostat=ios)
        if (ios /= 0) then
            print *, "Error in read_namelist:"
            print *, "iostat = ", ios
            error stop
        end if
        close (unit)

        call check_if_valid_namelist()

    end subroutine read_namelist

    subroutine check_if_valid_namelist()
        logical :: is_valid

        is_valid = .true.

        if (len(trim(log_file)) == 0) then
            print *, "field_file is empty!"
        end if
        if (len(trim(log_level)) == 0) then
            print *, "field_file is empty!"
        end if

        if (.not. is_valid) error stop

    end subroutine check_if_valid_namelist

end module logger_config
