program test_bmod_spectral
    use constants, only: dp, pi
    use bmod_spectral, only: bmod_spectral_t
    use utils, only: not_same

    implicit none

    integer, parameter :: n_theta = 8
    integer, parameter :: n_phi = 10
    integer, parameter :: nfp = 5
    real(dp), parameter :: reltol = 1e-12_dp
    real(dp), parameter :: abstol = 1e-12_dp

    type(bmod_spectral_t) :: spectrum
    real(dp) :: values(n_theta, n_phi)
    real(dp) :: theta, phi, found, expected
    integer :: j, k

    do k = 1, n_phi
        phi = 2.0_dp*pi*real(k - 1, dp)/real(nfp*n_phi, dp)
        do j = 1, n_theta
            theta = 2.0_dp*pi*real(j - 1, dp)/real(n_theta, dp)
            values(j, k) = analytic_bmod(theta, phi)
        end do
    end do

    call spectrum%build(values, nfp)

    call assert_value(spectrum, 0.17_dp, 0.03_dp)
    call assert_value(spectrum, 1.31_dp, 0.29_dp)
    call assert_value(spectrum, 5.01_dp, 0.88_dp)
    call assert_value(spectrum, 2.0_dp*pi + 0.43_dp, 2.0_dp*pi/nfp + 0.11_dp)

    call spectrum%destroy()
    if (spectrum%ready) then
        print *, "test_bmod_spectral failed: destroy kept spectrum ready"
        error stop
    end if

contains

    real(dp) function analytic_bmod(theta, phi) result(bmod)
        real(dp), intent(in) :: theta, phi

        bmod = 3.0_dp &
               + 0.75_dp*cos(2.0_dp*theta - 3.0_dp*real(nfp, dp)*phi) &
               - 0.50_dp*sin(3.0_dp*theta + 2.0_dp*real(nfp, dp)*phi) &
               + 0.20_dp*cos(4.0_dp*real(nfp, dp)*phi)
    end function analytic_bmod

    subroutine assert_value(spectrum, theta, phi)
        type(bmod_spectral_t), intent(in) :: spectrum
        real(dp), intent(in) :: theta, phi

        found = spectrum%evaluate(theta, phi)
        expected = analytic_bmod(theta, phi)
        if (not_same(found, expected, reltol_in=reltol, abstol_in=abstol)) then
            print *, "test_bmod_spectral failed"
            print *, "theta: ", theta
            print *, "phi: ", phi
            print *, "found: ", found
            print *, "expected: ", expected
            error stop
        end if
    end subroutine assert_value

end program test_bmod_spectral
