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
