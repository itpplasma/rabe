program test_integrate_wrapper
    use constants, only: dp
    use integrate, only: integrate_1d
    use utils, only: is_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-8_dp, abstol = 0.0_dp
    real(dp), dimension(2), parameter :: interval = (/0.0_dp, 1.0_dp/)
    real(dp), parameter :: constant_1 = 0.0_dp, constant_2 = 1.0_dp
    real(dp), parameter :: integral_1 = 1.0_dp, integral_2 = 2.0_dp

    real(dp) :: wrapper_constant
    real(dp) :: found_integral_1, found_integral_2

    wrapper_constant = constant_1
    call integrate_1d(polynome_wrapper, interval(1), interval(2), found_integral_1)

    if (is_same(integral_1, found_integral_1, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_wrapper failed: integral"
        print *, "found: ", found_integral_1
        print *, "analytic: ", integral_1
        error stop
    end if

    wrapper_constant = constant_2
    call integrate_1d(polynome_wrapper, interval(1), interval(2), found_integral_2)

    if (is_same(integral_2, found_integral_2, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_wrapper failed: integral"
        print *, "found: ", found_integral_2
        print *, "analytic: ", integral_2
        error stop
    end if

contains

    function polynome_wrapper(x)
        real(dp), intent(in) :: x
        real(dp) :: polynome_wrapper

        polynome_wrapper = polynome(x, wrapper_constant)
    end function polynome_wrapper

    function polynome(x, constant)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: constant
        real(dp) :: polynome

        polynome = 4*x**3 + constant
    end function polynome

end program test_integrate_wrapper
