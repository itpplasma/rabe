module error_handling
    implicit none

    public :: error_stop_unless_unsafe, set_unsafe_mode

    logical, private :: unsafe_mode = .false.

contains

    subroutine set_unsafe_mode(value)
        logical, intent(in) :: value
        unsafe_mode = value
    end subroutine set_unsafe_mode

    subroutine error_stop_unless_unsafe()
        if (unsafe_mode) return
        error stop
    end subroutine error_stop_unless_unsafe

end module error_handling
