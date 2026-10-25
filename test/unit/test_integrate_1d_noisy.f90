program test_integrate_1d_noisy
    ! Negative control for the roundoff acceptance in integrate_1d: a nearly
    ! cancelling integrand with deterministic noise of amplitude 1e-3 makes dqag
    ! return error flag 2 with an error above 1e-6*int |f|, so integrate_1d
    ! must still stop (registered with WILL_FAIL).
    use constants, only: dp, pi
    use integrate, only: integrate_1d

    implicit none

    real(dp) :: found_integral

    call integrate_1d(noisy, 0.0_dp, 2.0_dp*pi, found_integral)
    print *, "test_integrate_1d_noisy: no stop, result ", found_integral

contains

    function noisy(x)
        real(dp), intent(in) :: x
        real(dp) :: noisy, hash

        hash = sin(12345.678_dp*x)*43758.5453_dp
        noisy = sin(x) + 1.0e-3_dp*(hash - floor(hash) - 0.5_dp)
    end function noisy

end program test_integrate_1d_noisy
