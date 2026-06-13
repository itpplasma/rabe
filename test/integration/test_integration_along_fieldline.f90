program test_integration_along_fieldline
    use constants, only: dp, pi
    use utils, only: not_same
    use mock_field, only: mock_field_t
    use fieldline_mod, only: fieldline_t
    use field_checks, only: suspect_omnigenous_origin_not_minimum
    use make_fieldline, only: find_maxima_along_fieldline
    use make_fieldline, only: maxima_t
    use integrate, only: integrate_1d

    implicit none

    real(dp), parameter :: error_limit = 1e-7
    real(dp) :: reltol

    real(dp), parameter :: M_pol = 1.0_dp, N_tor = -4.0_dp
    real(dp), parameter :: nfp = max(1.0_dp, abs(N_tor))
    real(dp), parameter :: B_0 = 1.0_dp, B_amplitude = -0.5_dp
    real(dp), parameter :: integral = 2.0_dp*pi*B_0
    type(mock_field_t) :: field

    real(dp), parameter :: theta_0 = -3.0_dp*pi
    real(dp), parameter :: iota = -3.0_dp
    type(fieldline_t), dimension(1) :: fieldline

    type(maxima_t) :: maxima
    real(dp) :: interval(2)
    real(dp) :: found_integral

    call field%mock_field_init(M_pol, N_tor, B_0, B_amplitude)

    fieldline(1)%xi_0 = 0.0_dp
    fieldline(1)%theta_0 = theta_0
    fieldline(1)%iota = iota

    if (suspect_omnigenous_origin_not_minimum(field, M_pol, N_tor)) then
        print *, "error: The origin of the IDEAL omnigenous configuration"
        print *, "(theta=phi=0) must be a global and local minimum!"
        print *, "Origin of provided field suggests that this is not the case!"
        error stop
    end if
    fieldline%theta_0 = N_tor*fieldline%xi_0/nfp
    fieldline%phi_0 = M_pol*fieldline%xi_0/nfp

    interval = [0.0_dp, 4*pi] + fieldline(1)%phi_0
    call find_maxima_along_fieldline(field, fieldline(1), interval, maxima)
    if (maxima%n == 2) then
        fieldline(1)%phi_max = maxima%phi(1:2)
    else
        print *, "-------------------------------------------------------------"
        print *, "test_integration_along_fieldline failed: maxima"
        print *, "Found ", maxima%n, " maxima, expected 2!"
        print *, "maxima: ", maxima%phi(1:max(1, maxima%n))
        error stop
    end if
    if (maxval(maxima%phi_error) > error_limit) then
        print *, "-------------------------------------------------------------"
        print *, "test_integration_along_fieldline failed: maxima"
        print *, "Achieved error is too large: ", maxima%phi_error
        print *, "Error limit: ", error_limit
        error stop
    end if

    call integrate_1d(B_mod_along_fieldline, &
                      fieldline(1)%phi_max(1), &
                      fieldline(1)%phi_max(2), &
                      found_integral)

    reltol = maxval(maxima%phi_error)*2.0_dp
    if (not_same(integral, found_integral, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integration_along_fieldline failed: integral"
        print *, "found: ", found_integral
        print *, "expected: ", integral
        error stop
    end if

contains
    function B_mod_along_fieldline(phi) result(B_mod)
        real(dp), intent(in) :: phi
        real(dp) :: B_mod

        real(dp) :: theta

        theta = (phi - fieldline(1)%phi_0)*fieldline(1)%iota + fieldline(1)%theta_0
        call field%compute_B_mod(theta, phi, B_mod)
    end function B_mod_along_fieldline

end program test_integration_along_fieldline
