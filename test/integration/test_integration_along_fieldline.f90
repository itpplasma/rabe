program test_integration_along_fieldline
    use constants, only: dp, pi
    use utils, only: is_same
    use mock_field, only: mock_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_mod, only: set_fieldline_phi_0_to_mode_minimum
    use fieldline_mod, only: find_maxima_along_fieldline

    implicit none

    real(dp), parameter :: reltol = 1e-2

    real(dp), parameter :: theta_mode = 1.0_dp, phi_mode = -4.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, B_amplitude = 0.5_dp
    type(mock_field_t) :: field

    real(dp), parameter :: theta_0 = -3.0_dp*pi
    real(dp), parameter :: iota = -3.0_dp
    integer, parameter :: n_maxima = 2
    type(fieldline_t), dimension(1) :: fieldline

    real(dp) :: interval(2)
    real(dp) :: integral, found_integral

    call field%mock_field_init(theta_mode, phi_mode, B_0, B_amplitude)

    fieldline(1)%theta_0 = theta_0
    fieldline(1)%iota = iota
    allocate (fieldline(1)%phi_max(n_maxima))

    call set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, fieldline)

    interval = (/0.0_dp, 4*pi/) + fieldline(1)%phi_0
    call find_maxima_along_fieldline(field, fieldline(1), interval, &
                                     fieldline(1)%phi_max)

    found_integral = 0.0_dp

    integral = 2.0_dp*pi*B_0
    if (is_same(integral, found_integral, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integration_along_fieldline failed: integral"
        print *, "found: ", found_integral
        print *, "expected: ", integral
        error stop
    end if

end program test_integration_along_fieldline
