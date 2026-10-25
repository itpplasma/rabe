program test_integrate_1d_cancellation
    ! A nearly cancelling integral, int_0^{2 pi} (sin(x) + delta) dx = 2 pi delta,
    ! cannot reach the relative tolerance 1e-6 on |result| because quadpack's
    ! error estimate stays at the roundoff of int |f| = 4. dqag reports error
    ! flag 2, which integrate_1d must accept since the error is below 1e-6*4.
    use constants, only: dp, pi
    use integrate, only: integrate_1d

    implicit none

    real(dp), parameter :: delta = 1.0e-12_dp
    real(dp), parameter :: abstol = 1.0e-13_dp
    real(dp) :: found_integral

    call integrate_1d(nearly_cancelling, 0.0_dp, 2.0_dp*pi, found_integral)

    if (abs(found_integral - 2.0_dp*pi*delta) > abstol) then
        print *, "test_integrate_1d_cancellation failed"
        print *, "found: ", found_integral
        print *, "analytic: ", 2.0_dp*pi*delta
        error stop
    end if

contains

    function nearly_cancelling(x)
        real(dp), intent(in) :: x
        real(dp) :: nearly_cancelling

        nearly_cancelling = sin(x) + delta
    end function nearly_cancelling

end program test_integrate_1d_cancellation
