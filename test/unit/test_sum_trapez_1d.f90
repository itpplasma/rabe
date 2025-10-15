program test_sum_trapez_1d
    use constants, only: dp, pi
    use integrate, only: sum_trapez_1d

    implicit none

    real(dp), dimension(2), parameter :: interval_0 = [-1.0_dp, 1.0_dp]
    real(dp), parameter :: integral_0 = 2.0_dp
    real(dp), dimension(2), parameter :: interval_1 = [0.0_dp, 0.5_dp*pi]
    real(dp), parameter :: integral_1 = 2.0_dp*(sqrt(2.0_dp) - 1.0_dp)
    real(dp), dimension(2), parameter :: interval_2 = [-pi, pi]
    real(dp), parameter :: integral_2 = 4.0_dp*sqrt(2.0_dp)
    real(dp), dimension(2), parameter :: interval_3 = [-pi, pi]
    real(dp), parameter :: eps = 0.5_dp
    real(dp), parameter :: integral_3 = integral_2*(1.0_dp + eps/3.0_dp)

    logical :: failed_test

    failed_test = .false.

    print *, "integral 0"
    if (err_not_decreasing_with_refinement(trial_func_0, interval_0, integral_0)) then
        failed_test = .true.
    end if
    print *, "integral 1"
    if (err_not_decreasing_with_refinement(trial_func_1, interval_1, integral_1)) then
        failed_test = .true.
    end if
    print *, "integral 2"
    if (err_not_decreasing_with_refinement(trial_func_2, interval_2, integral_2)) then
        failed_test = .true.
    end if
    print *, "integral 3"
    if (err_not_decreasing_with_refinement(trial_func_3, interval_3, integral_3)) then
        failed_test = .true.
    end if

    if (failed_test) error stop

contains

    function err_not_decreasing_with_refinement(trial_func, interval, integral)
        use utils, only: not_same, linspace
        interface
            function func(x)
                import :: dp
                real(dp), intent(in) :: x
                real(dp) :: func
            end function func
        end interface
        procedure(func) :: trial_func
        real(dp), intent(in) :: interval(2), integral
        logical :: err_not_decreasing_with_refinement

        integer, parameter :: n_x = 100, n_x_refined = 1000
        real(dp), dimension(:), allocatable :: x, y
        real(dp) :: err, err_expect
        integer :: current

        real(dp) :: found_integral

        err_not_decreasing_with_refinement = .false.

        allocate (x(n_x), y(n_x))
        call linspace(interval(1), interval(2), n_x, x)
        do current = 1, size(x)
            y(current) = trial_func(x(current))
        end do
        found_integral = sum_trapez_1d(x, y)
        deallocate (x, y)
        err = abs(found_integral - integral)
        allocate (x(n_x_refined), y(n_x_refined))
        call linspace(interval(1), interval(2), n_x_refined, x)
        do current = 1, size(x)
            y(current) = trial_func(x(current))
        end do
        found_integral = sum_trapez_1d(x, y)
        deallocate (x, y)
        err_expect = err*(real(n_x, kind=dp)/real(n_x_refined, kind=dp))**2.0_dp

        if (not_same(integral, found_integral, abstol_in=err_expect)) then
            print *, "-------------------------------------------------------------"
            print *, "test_sum_trapez_1d failed:"
            print *, "found: ", found_integral
            print *, "analytic: ", integral
            print *, "error: ", abs(found_integral - integral)
            print *, "expected error: ", err_expect
            err_not_decreasing_with_refinement = .true.
        end if
    end function err_not_decreasing_with_refinement

    function trial_func_0(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_0

        trial_func_0 = 3.0_dp*x**2
    end function trial_func_0

    function trial_func_1(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_1

        trial_func_1 = sqrt(1.0_dp - sin(x))
    end function trial_func_1

    function trial_func_2(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_2

        trial_func_2 = sqrt(1.0_dp + cos(x))
    end function trial_func_2

    function trial_func_3(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_3

        trial_func_3 = trial_func_2(x)*(1.0_dp + eps*cos(x))
    end function trial_func_3

end program test_sum_trapez_1d
