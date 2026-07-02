module error_handling
    implicit none

    private
    public :: unsafe_mode, read_error_handling_config
    public :: failed_sanity_check, set_unsafe_mode
    public :: reset_failed_check_counter, did_fail_any_sanity_check

    logical, protected :: unsafe_mode = .false.
    integer, private :: failed_check_counter = 0

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

    !>
    !! \brief Toggle unsafe mode for sanity checks.
    !!
    !! \details When .true., a failed check increments the counter instead of halting,
    !! allowing surface-by-surface iteration to continue. Call
    !! did_fail_any_sanity_check to detect failures.
    !<
    subroutine set_unsafe_mode(value)
        logical, intent(in) :: value
        unsafe_mode = value
    end subroutine set_unsafe_mode

    subroutine failed_sanity_check()
        if (unsafe_mode) then
            failed_check_counter = failed_check_counter + 1
            return
        else
            error stop
        end if
    end subroutine failed_sanity_check

    !>
    !! \brief Reset the failed sanity check counter to zero.
    !<
    subroutine reset_failed_check_counter()
        failed_check_counter = 0
    end subroutine reset_failed_check_counter

    !>
    !! \brief Return .true. if any sanity check has failed since the last reset.
    !<
    function did_fail_any_sanity_check()
        logical :: did_fail_any_sanity_check
        did_fail_any_sanity_check = (failed_check_counter > 0)
    end function did_fail_any_sanity_check

end module error_handling
