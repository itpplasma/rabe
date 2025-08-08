program test_fieldline_integral
    use constants, only: dp, pi
    use utils, only: not_same
    use mock_field, only: mock_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_mod, only: set_fieldline_phi_0_to_mode_minimum
    use fieldline_mod, only: find_maxima_along_fieldline
    use integrate, only: integrate_1d

    implicit none

    real(dp), parameter :: reltol = 1e-3 !retol >= retol of find maxima
    integer, parameter :: n_steps = 10000 !retol ~ interval/n_steps

    real(dp), parameter :: theta_mode = 1.0_dp, phi_mode = -4.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, B_amplitude = 0.5_dp
    real(dp), parameter :: integral = 2.0_dp*pi*B_0
    type(mock_field_t) :: field

    real(dp), parameter :: theta_0 = -3.0_dp*pi
    real(dp), parameter :: iota = -3.0_dp
    type(fieldline_t), dimension(1) :: fieldline

    real(dp) :: interval(2)
    real(dp) :: found_integral

    call field%mock_field_init(theta_mode, phi_mode, B_0, B_amplitude)

    fieldline(1)%theta_0 = theta_0
    fieldline(1)%iota = iota

    call set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, fieldline)

    interval = [0.0_dp, 4.0_dp*pi] + fieldline(1)%phi_0
    call find_maxima_along_fieldline(field, fieldline(1), interval, n_steps)

    call integrate_1d(B_mod_along_fieldline, &
                      fieldline(1)%phi_max(1), &
                      fieldline(1)%phi_max(2), &
                      found_integral)

    if (not_same(integral, found_integral, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_fieldline_integral failed: integral"
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

end program test_fieldline_integral
