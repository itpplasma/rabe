program test_calc_periodic_dydx
    use constants, only: dp, pi
    use utils, only: not_same, linspace
    use shaing_callen_remainder, only: calc_periodic_dydx
    implicit none

    integer, parameter :: ns(5) = [50, 100, 200, 400, 800]
    real(dp), parameter :: reltol = 0.0_dp
    real(dp), parameter :: range = 2.0_dp*pi
    real(dp) :: abstol

    real(dp), dimension(:), allocatable :: x, y, dydx, temp
    integer :: this
    integer :: n

    logical :: test_failed

    test_failed = .false.
    do this = 1, size(ns)
        n = ns(this)
        abstol = (range)**2.0_dp/6.0_dp/real(n, kind=dp)**2.0_dp
        allocate (x(n), y(n), dydx(n), temp(n + 1))
        call linspace(-pi, pi, n + 1, temp)
        x = temp(1:n)
        y = trial_function(x)
        dydx = calc_periodic_dydx(x, y)
        if (not_same(dydx, trial_derivative(x), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_periodic_dydx failed:"
            print *, "found: ", dydx
            print *, "analytic: ", trial_derivative(x)
            print *, "max abs error: ", maxval(abs(dydx - trial_derivative(x)))
            print *, "for number of points = ", n
            print *, "expected error: ", abstol
            test_failed = .true.
        end if
        deallocate (x, y, dydx, temp)
    end do

    if (test_failed) error stop

contains
    elemental function trial_function(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp), parameter :: phase = 0.2_dp*pi

        y = sin(x + phase)
    end function trial_function

    elemental function trial_derivative(x) result(dydx)
        real(dp), intent(in) :: x
        real(dp) :: dydx
        real(dp), parameter :: phase = 0.2_dp*pi

        dydx = cos(x + phase)
    end function trial_derivative

end program test_calc_periodic_dydx
