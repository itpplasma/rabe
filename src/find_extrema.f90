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

    subroutine find_local_minima(func, interval, location, n_steps_in)
        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        real(dp), dimension(:), intent(out) :: location
        integer, intent(in), optional :: n_steps_in

        integer :: n_steps

        if (present(n_steps_in)) then
            n_steps = n_steps_in
        else
            n_steps = 1000
        end if

        call find_local_maxima(negative_func, interval, location, n_steps)

    contains
        subroutine negative_func(x, value)
            use constants, only: dp
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(:), intent(out) :: value
            call func(x, value)
            value = -value
        end subroutine negative_func

    end subroutine find_local_minima

    subroutine find_local_maxima(func, interval, location, n_steps_in)
        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        real(dp), dimension(:), intent(out) :: location
        integer, intent(in), optional :: n_steps_in

        integer :: n_maxima, n_steps
        real(dp), dimension(:), allocatable :: x, value
        integer :: current_maximum, current_location

        if (present(n_steps_in)) then
            n_steps = n_steps_in
        else
            n_steps = 1000
        end if
        allocate (x(n_steps), value(n_steps))
        call linspace(interval(1), interval(2), n_steps, x)
        call func(x, value)

        n_maxima = size(location, dim=1)
        current_maximum = 0
        do current_location = 2, n_steps - 1
            if (value(current_location - 1) < value(current_location)) then
                if (value(current_location + 1) < value(current_location)) then
                    current_maximum = current_maximum + 1
                    location(current_maximum) = x(current_location)
                    if (current_maximum .eq. n_maxima) exit
                end if
            end if
        end do

        if (current_maximum < n_maxima) then
            print *, "find_local_maxima: found less maxima then expected"
            print *, "found: ", current_maximum
            print *, "expected: ", n_maxima
            error stop
        end if

    end subroutine find_local_maxima

    function find_global_extrema(func, interval, n_steps_in) result(extrema)
        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        integer, intent(in), optional :: n_steps_in
        real(dp) :: extrema(2)

        integer :: n_steps
        real(dp), dimension(:), allocatable :: x, value

        if (present(n_steps_in)) then
            n_steps = n_steps_in
        else
            n_steps = 1000
        end if

        allocate (x(n_steps), value(n_steps))
        call linspace(interval(1), interval(2), n_steps, x)
        call func(x, value)
        extrema(1) = minval(value)
        extrema(2) = maxval(value)
    end function find_global_extrema

    subroutine linspace(a, b, n, x)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        real(dp), dimension(:), intent(out) :: x

        real(dp) :: dx
        integer :: i

        dx = (b - a)/(n - 1)
        do i = 1, n
            x(i) = a + (i - 1)*dx
        end do
    end subroutine linspace

end module find_extrema
