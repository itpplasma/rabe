program test_find_extrema
use, intrinsic :: iso_fortran_env, only: dp => real64
use find_extrema, only: find_global_extrema

implicit none

real(dp) :: pi = 3.14159
real(dp) :: tol = 1e-3

call test_find_global_extrema()


contains


subroutine test_find_global_extrema
    real(dp) :: extrema(2), expected_extrema(2)
    real(dp) :: interval(2)

    interval = [0.0_dp, 2.0_dp*pi]
    expected_extrema = (/-1.0_dp, 1.0_dp/)
    extrema = find_global_extrema(test_func, interval)

    if(any(abs(extrema/expected_extrema - 1) > tol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_find_global_extrema"
        error stop
    endif
end subroutine test_find_global_extrema

subroutine test_func(x, value)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(out) :: value

    value = sin(x)
end subroutine test_func

end program test_find_extrema