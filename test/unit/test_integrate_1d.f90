program test_integrate_1d
    use constants, only: dp
    use integrate, only: integrate_1d
    use utils, only: is_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-8_dp, abstol = 0.0_dp
    real(dp), dimension(2), parameter :: interval = (/0.0_dp, 1.0_dp/)
    real(dp), parameter :: integral = 1.0_dp

    real(dp) :: found_integral

    call integrate_1d(polynome, interval(1), interval(2), found_integral)

    if (is_same(integral, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d failed: integral"
        print *, "found: ", found_integral
        print *, "analytic: ", integral
        error stop
    end if

contains

    function polynome(x)
        real(dp), intent(in) :: x
        real(dp) :: polynome

        polynome = 4*x**3
    end function polynome

end program test_integrate_1d
