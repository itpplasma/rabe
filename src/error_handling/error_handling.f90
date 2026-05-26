module error_handling
    implicit none

    public :: error_stop_unless_unsafe, read_error_handling_config

    logical, protected :: unsafe_mode

contains

    subroutine read_error_handling_config(config_file)
        character(len=*), intent(in) :: config_file
        namelist /error_handling/ unsafe_mode
        integer :: ios, unit

        unsafe_mode = .false.
        open (newunit=unit, file=config_file, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print *, "read_error_handling_config: ", &
                "failed to open error handling config file: ", trim(config_file)
            error stop
        end if

        read (unit, nml=error_handling, iostat=ios)
        if (ios > 0) then
            print *, "read_error_handling_config: error reading namelist, iostat: ", ios
            error stop
        end if
        close (unit)

    end subroutine read_error_handling_config

    subroutine error_stop_unless_unsafe()
        if (unsafe_mode) return
        error stop
    end subroutine error_stop_unless_unsafe

end module error_handling
