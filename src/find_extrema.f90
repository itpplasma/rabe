module find_extrema
    use constants, only: dp
    implicit none

    abstract interface
        subroutine func1d(x, value)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(:), intent(out) :: value
        end subroutine func1d
    end interface

contains

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
            n_steps = 100
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
