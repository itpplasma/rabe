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
        end if

        error_limit = abs(result_quadpack)*rel_error_tol_quadpack &
                      + abs_error_tol_quadpack

        if (abs_error > error_limit) then
            print *, "Integration warning: integration error =", abs_error
            print *, "bigger than required ", error_limit
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
