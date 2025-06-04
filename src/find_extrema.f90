module find_extrema
    use constants, only: dp
    implicit none

    abstract interface
        subroutine func1d(x, value)
            use constants, only: dp
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(:), intent(out) :: value
        end subroutine func1d
    end interface

contains

    subroutine find_local_minima(func, interval, location, abstol)
        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        real(dp), dimension(:), intent(out) :: location
        real(dp), intent(in), optional :: abstol

        call find_local_maxima(negative_func, interval, location, abstol)

    contains
        subroutine negative_func(x, value)
            use constants, only: dp
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(:), intent(out) :: value
            call func(x, value)
            value = -value
        end subroutine negative_func

    end subroutine find_local_minima

    recursive subroutine find_local_maxima(func, interval, location, abstol)
        use utils, only: linspace

        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        real(dp), dimension(:), intent(out) :: location
        real(dp), intent(in), optional :: abstol

        integer :: n_maxima, n_steps
        integer, parameter :: n_max = 10000
        real(dp), dimension(:), allocatable :: x, value
        integer :: current_maximum, current_location

        logical :: do_recursion
        real(dp) :: error
        real(dp) :: subinterval(2)
        real(dp), dimension(1) :: sublocation

        if (present(abstol)) then
            n_steps = int(abs(interval(2) - interval(1))/abstol) + 2
        else
            n_steps = 1000
        end if

        if (n_steps <= 0 .or. n_steps >= n_max) then
            n_steps = 1000
            do_recursion = .true.
        else
            do_recursion = .false.
        end if

        allocate (x(n_steps), value(n_steps))
        call linspace(interval(1), interval(2), n_steps, x)
        call func(x, value)

        error = abs(x(2) - x(1))

        n_maxima = size(location, dim=1)
        current_maximum = 0
        do current_location = 2, n_steps - 1
            if (value(current_location - 1) <= value(current_location)) then
                if (value(current_location + 1) <= value(current_location)) then
                    current_maximum = current_maximum + 1
                    if (do_recursion) then
                        if (error > abstol) then
                            subinterval(1) = x(current_location - 1)
                            subinterval(2) = x(current_location + 1)
                            call find_local_maxima(func, subinterval, sublocation, abstol)
                            location(current_maximum) = sublocation(1)
                        else
                            location(current_maximum) = x(current_location)
                        end if
                    else
                        location(current_maximum) = x(current_location)
                    end if
                    if (current_maximum .eq. n_maxima) exit
                end if
            end if
        end do

        if (current_maximum < n_maxima) then
            print *, "find_local_maxima: found less maxima then expected"
            print *, "found: ", current_maximum, "out of ", n_maxima
            print *, "location: ", location
            error stop
        end if

    end subroutine find_local_maxima

    function find_global_extrema(func, interval, abstol) result(extrema)
        use utils, only: linspace

        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        integer, intent(in), optional :: abstol
        real(dp) :: extrema(2)

        integer :: n_steps
        real(dp), dimension(:), allocatable :: x, value

        if (present(abstol)) then
            n_steps = int(abs(interval(2) - interval(1))/abstol) + 2
        else
            n_steps = 1000
        end if

        allocate (x(n_steps), value(n_steps))
        call linspace(interval(1), interval(2), n_steps, x)
        call func(x, value)
        extrema(1) = minval(value)
        extrema(2) = maxval(value)
    end function find_global_extrema

end module find_extrema
