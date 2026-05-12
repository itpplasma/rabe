module integrate
    use constants, only: dp
    use odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none

    abstract interface
        function integrand(x)
            use constants, only: dp
            real(dp), intent(in) :: x
            real(dp) :: integrand
        end function integrand
    end interface

contains

    subroutine integrate_1d(f, a, b, result)
        procedure(integrand) :: f
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: result

        real(dp), parameter :: eps = 1.0e-10_dp
        real(dp), parameter :: eps_rough = 1.0e-3_dp
        real(dp) :: y(1), y_offset

        ! Round 1: cheap pass to estimate the integral scale.
        y(1) = eps_rough
        call odeint_allroutines(y, 1, a, b, eps_rough, integral_derivs)
        y_offset = max(abs(y(1) - eps_rough), eps)

        ! Round 2: accurate pass with y_offset calibrated to the integral magnitude,
        ! so y_scale ~ |integral| throughout and relative error is ~eps.
        y(1) = y_offset
        call odeint_allroutines(y, 1, a, b, eps, integral_derivs)

        if (odeint_has_failed()) then
            print *, "Integration failed!"
            error stop
        end if

        result = y(1) - y_offset

        if (ieee_is_nan(result)) then
            print *, "Integration result is NaN!"
            error stop
        end if

    contains

        subroutine integral_derivs(x, yy, dydx)
            real(dp), intent(in) :: x
            real(dp), intent(in) :: yy(:)
            real(dp), intent(out) :: dydx(:)
            dydx(1) = f(x)
        end subroutine integral_derivs

    end subroutine integrate_1d

    subroutine integrate_1d_substituted(f, a, b, result)
        procedure(integrand) :: f
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: result

        real(dp) :: sub_a, sub_b
        real(dp) :: left_result, right_result

        if (b < a) then
            print *, "Error in integrate_1d_substituted: b < a"
            error stop
        end if

        ! integration by globally adaptive interval subdivision (quadpack import)
        ! qags struggels with functions of form sqrt(x-a) and sqrt(b-x)
        ! Terefore we split the integral in the middle and substitute
        ! - left integral: t**2 = x - a
        ! - right integral: t**2 = b - x
        ! resulting in new limits
        ! - left integral: 0 to sqrt((b - a)*0.5_dp)
        ! - right integral: 0 to sqrt((b - a)*0.5_dp) (after absorbing [-] from trafo)
        sub_a = 0.0_dp
        sub_b = sqrt((b - a)*0.5_dp)
        call integrate_1d(left_substituted_integrand, sub_a, sub_b, left_result)
        call integrate_1d(right_substituted_integrand, sub_a, sub_b, right_result)

        result = left_result + right_result

    contains
        function left_substituted_integrand(t)
            real(dp), intent(in) :: t
            real(dp) :: left_substituted_integrand

            real(dp) :: x

            x = t**2.0_dp + a
            left_substituted_integrand = f(x)*2.0_dp*t
        end function left_substituted_integrand

        function right_substituted_integrand(t)
            real(dp), intent(in) :: t
            real(dp) :: right_substituted_integrand

            real(dp) :: x

            x = b - t**2.0_dp
            right_substituted_integrand = f(x)*2.0_dp*t
        end function right_substituted_integrand

    end subroutine integrate_1d_substituted

    function sum_trapez_1d(x, y)
        real(dp), dimension(:), intent(in) :: x, y
        real(dp) :: sum_trapez_1d

        integer :: n

        n = size(x)
        if (.not. n == size(y)) then
            print *, "Error in sum_trapez_1d: x and y need to have same length!"
            error stop
        end if
        if (n < 2) then
            print *, "Error in sum_trapez_1d: len(x) must be at least 2!"
            error stop
        end if

        sum_trapez_1d = 0.5_dp*sum((x(2:n) - x(1:n - 1))*(y(2:n) + y(1:n - 1)))
    end function sum_trapez_1d

end module integrate
