module error_handling
    implicit none

    public :: failed_sanity_check, set_unsafe_mode
    public :: did_fail_any_sanity_check

    logical, private :: unsafe_mode = .false.
    integer, private :: failed_check_counter = 0

contains

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

    subroutine reset_failed_check_counter()
        failed_check_counter = 0
    end subroutine reset_failed_check_counter

    function did_fail_any_sanity_check()
        logical :: did_fail_any_sanity_check
        did_fail_any_sanity_check = (failed_check_counter > 0)
    end function did_fail_any_sanity_check

end module error_handling
