program test_integrate_1d_substituted
    use constants, only: dp, pi
    use integrate, only: integrate_1d_substituted, integrate_1d
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-8_dp, abstol = 0.0_dp
    real(dp), dimension(2), parameter :: interval = (/-1.0_dp, 1.0_dp/)
    real(dp), parameter :: integral = 0.5_dp*pi

    real(dp) :: found_integral

    call integrate_1d_substituted(trial_func, interval(1), interval(2), found_integral)

    if (not_same(integral, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d failed: integral"
        print *, "found: ", found_integral
        print *, "analytic: ", integral
        error stop
    end if

contains

    function trial_func(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func

        trial_func = sqrt(1 - x**2)
    end function trial_func

end program test_integrate_1d_substituted
