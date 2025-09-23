module integrate
    use constants, only: dp

    implicit none

    integer, parameter :: quadpack = 8

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

        real(quadpack), parameter :: abs_error_tol_quadpack = 0.0_dp
        real(quadpack), parameter :: rel_error_tol_quadpack = 1.0e-6
        integer, parameter :: order_key = 6
        real(quadpack) :: a_quadpack, b_quadpack, result_quadpack
        real(quadpack) :: abs_error
        integer :: error_flag
        real(quadpack) :: error_limit

        integer :: neval_dummy

        ! quadpack operates in and requires real(8)
        a_quadpack = convert_to_quadpack(a)
        b_quadpack = convert_to_quadpack(b)

        ! integration by globally adaptive interval subdivision (quadpack import)
        call qag(quadpack_integrand, &
                 a_quadpack, &
                 b_quadpack, &
                 abs_error_tol_quadpack, &
                 rel_error_tol_quadpack, &
                 order_key, &
                 result_quadpack, &
                 abs_error, &
                 neval_dummy, &
                 error_flag)

        error_limit = abs(result_quadpack)*rel_error_tol_quadpack &
                      + abs_error_tol_quadpack

        if (abs_error > error_limit) then
            print *, "Integration warning: absolute error =", abs_error
            print *, "bigger than required ", error_limit
            print *, "relative error", abs_error/abs(result_quadpack)
            error stop
        end if

        if (error_flag /= 0) then
            print *, "Integration warning: error =", error_flag
            error stop
        end if

        result = convert_to_dp(result_quadpack)

    contains

        ! the function input also needs to be real(8) for quadpack
        function quadpack_integrand(x_quadpack)
            real(quadpack), intent(in) :: x_quadpack
            real(quadpack) :: quadpack_integrand

            real(dp) :: x_dp

            x_dp = real(x_quadpack, kind=dp)
            quadpack_integrand = real(f(x_dp), kind=quadpack)
        end function quadpack_integrand

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

    function convert_to_quadpack(val_dp) result(val_quadpack)
        real(dp), intent(in) :: val_dp
        real(quadpack) :: val_quadpack

        val_quadpack = real(val_dp, kind=quadpack)
    end function convert_to_quadpack

    function convert_to_dp(val_quadpack) result(val_dp)
        real(quadpack), intent(in) :: val_quadpack
        real(dp) :: val_dp

        val_dp = real(val_quadpack, kind=dp)
    end function convert_to_dp

end module integrate
