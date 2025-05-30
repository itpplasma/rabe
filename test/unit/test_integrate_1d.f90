program test_integrate_1d
    use constants, only: dp, pi
    use integrate, only: integrate_1d
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-8_dp, abstol = 0.0_dp
    real(dp), dimension(2), parameter :: interval_polynom = (/0.0_dp, 1.0_dp/)
    real(dp), parameter :: integral_polynom = 1.0_dp
    real(dp), dimension(2), parameter :: interval_not_polynom = (/0.0_dp, 0.5_dp*pi/)
    real(dp), parameter :: integral_not_polynom = 10.0_dp

    real(dp) :: found_integral

    call integrate_1d(polynome, &
                      interval_polynom(1), &
                      interval_polynom(2), &
                      found_integral)

    if (not_same(integral_polynom, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d failed: integral of polynom"
        print *, "found: ", found_integral
        print *, "analytic: ", integral_polynom
        error stop
    end if

    call integrate_1d(not_polynome, &
                      interval_not_polynom(1), &
                      interval_not_polynom(2), &
                      found_integral)

    if (not_same(integral_not_polynom, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d failed: integral of not-polynom"
        print *, "found: ", found_integral
        print *, "analytic: ", integral_not_polynom
        error stop
    end if

contains

    function polynome(x)
        real(dp), intent(in) :: x
        real(dp) :: polynome

        polynome = 4*x**3
    end function polynome

    function not_polynome(x)
        real(dp), intent(in) :: x
        real(dp) :: not_polynome

        not_polynome = 10.0_dp*sin(x)
    end function not_polynome

end program test_integrate_1d
