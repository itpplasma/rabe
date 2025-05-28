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

        real(quadpack), parameter :: abs_error_tol_quadpack = 1.0e-8
        real(quadpack), parameter :: rel_error_tol_quadpack = 1.0e-8
        real(quadpack) :: a_quadpack, b_quadpack, result_quadpack
        real(quadpack) :: abs_error
        integer :: error_flag
        real(quadpack) :: error_limit

        integer :: neval_dummy

        ! quadpack operates in and requires real(8)
        a_quadpack = convert_to_quadpack(a)
        b_quadpack = convert_to_quadpack(b)

        ! integration by globally adaptive interval subdivision (quadpack import)
        call qags(quadpack_integrand, &
                  a_quadpack, &
                  b_quadpack, &
                  abs_error_tol_quadpack, &
                  rel_error_tol_quadpack, &
                  result_quadpack, &
                  abs_error, &
                  neval_dummy, &
                  error_flag)

        if (error_flag /= 0) then
            print *, "Integration warning: error =", error_flag
            error stop
        end if

        error_limit = abs(result_quadpack)*rel_error_tol_quadpack &
                      + abs_error_tol_quadpack

        if (abs_error > error_limit) then
            print *, "Integration warning: integration error =", abs_error
            print *, "bigger than required ", error_limit
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

        real(dp) :: interval_middle

        real(quadpack), parameter :: abs_error_tol_quadpack = 1.0e-8
        real(quadpack), parameter :: rel_error_tol_quadpack = 1.0e-8
        real(quadpack) :: sub_a_quadpack, sub_b_quadpack
        real(quadpack) :: left_result_quadpack, right_result_quadpack
        real(quadpack) :: left_abs_error, right_abs_error
        integer :: left_error_flag, right_error_flag
        real(quadpack) :: left_error_limit, right_error_limit

        integer :: neval_dummy

        if (b < a) then
            print *, "Error in integrate_1d_substituted: b < a"
            error stop
        end if
        interval_middle = (b - a)*0.5_dp

        ! quadpack operates in and requires real(8)
        sub_a_quadpack = convert_to_quadpack(0.0_dp)
        sub_b_quadpack = convert_to_quadpack(sqrt(interval_middle))

        ! integration by globally adaptive interval subdivision (quadpack import)
        ! qags struggels with functions of form sqrt(x-a) and sqrt(b-x)
        ! Terefore we split the integral in two parts and substitute
        ! - left integral: t**2 = x - a
        ! - right integral: t**2 = b - x
        call qags(left_substituted_integrand, &
                  sub_a_quadpack, &
                  sub_b_quadpack, &
                  abs_error_tol_quadpack, &
                  rel_error_tol_quadpack, &
                  left_result_quadpack, &
                  left_abs_error, &
                  neval_dummy, &
                  left_error_flag)

        call qags(right_substituted_integrand, &
                  sub_a_quadpack, &
                  sub_b_quadpack, &
                  abs_error_tol_quadpack, &
                  rel_error_tol_quadpack, &
                  right_result_quadpack, &
                  right_abs_error, &
                  neval_dummy, &
                  right_error_flag)

        if (left_error_flag /= 0 .or. right_error_flag /= 0) then
            print *, "Integration warning: error =", left_error_flag
            print *, "Integration warning: error =", right_error_flag
            error stop
        end if

        left_error_limit = max(abs(left_result_quadpack)*rel_error_tol_quadpack, &
                               abs_error_tol_quadpack)
        right_error_limit = max(abs(right_result_quadpack)*rel_error_tol_quadpack, &
                                abs_error_tol_quadpack)

        if (left_abs_error > left_error_limit .or. &
            right_abs_error > right_error_limit) then
            print *, "Integration warning: too big integration error"
            print *, "reached in left: ", left_abs_error
            print *, "limit: ", left_error_limit
            print *, "reached in right: ", right_abs_error
            print *, "limit: ", right_error_limit
            error stop
        end if

        result = convert_to_dp(left_result_quadpack + right_result_quadpack)

    contains
        function left_substituted_integrand(t_quadpack)
            real(quadpack), intent(in) :: t_quadpack
            real(quadpack) :: left_substituted_integrand

            real(dp) :: t_dp, x_dp

            t_dp = real(t_quadpack, kind=dp)
            x_dp = t_dp**2 + a
            left_substituted_integrand = real(f(x_dp)*2.0_dp*t_dp, kind=quadpack)
        end function left_substituted_integrand

        function right_substituted_integrand(t_quadpack)
            real(quadpack), intent(in) :: t_quadpack
            real(quadpack) :: right_substituted_integrand

            real(dp) :: t_dp, x_dp

            t_dp = real(t_quadpack, kind=dp)
            x_dp = b - t_dp**2
            right_substituted_integrand = real(f(x_dp)*2.0_dp*t_dp, kind=quadpack)
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
