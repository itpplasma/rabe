program test_find_global_maximum
    use constants, only: dp, pi
    use utils, only: is_same, linspace
    use mock_field, only: mock_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_mod, only: set_fieldline_phi_0_to_mode_minimum
    use fieldline_mod, only: find_maxima_along_fieldline

    implicit none

    real(dp), parameter :: reltol = 1e-5 !retol >= retol of find maxima
    integer, parameter :: n_steps = 1000 !retol ~ interval/n_steps

    real(dp), parameter :: theta_mode = 1.0_dp, phi_mode = -4.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, B_amplitude = 0.5_dp
    real(dp), parameter :: global_B_max = B_0 + B_amplitude
    type(mock_field_t) :: field

    integer, parameter :: n_fieldlines = 10
    real(dp), parameter :: iota = -3.0_dp

    real(dp), dimension(n_fieldlines) :: theta_0
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp) :: interval(2)
    real(dp) :: found_global_B_max, fieldline_B_max
    integer :: current

    call field%mock_field_init(theta_mode, phi_mode, B_0, B_amplitude)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, theta_0)
    fieldlines(:)%theta_0 = theta_0
    fieldlines(:)%iota = iota

    call set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, fieldlines)

    found_global_B_max = 0.0_dp
    do current = 1, n_fieldlines
        interval = (/0.0_dp, 4*pi/) + fieldlines(current)%phi_0
        call find_maxima_along_fieldline(field, fieldlines(current), &
                                         interval, n_steps=n_steps)
        fieldline_B_max = max(fieldlines(current)%B_max(1), &
                              fieldlines(current)%B_max(2))
        if (fieldline_B_max .gt. found_global_B_max) then
            found_global_B_max = fieldline_B_max
        end if
    end do

    if (is_same(global_B_max, found_global_B_max, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_find_global_maximum failed: global_B_max"
        print *, "found: ", found_global_B_max
        print *, "expected: ", global_B_max
        error stop
    end if

end program test_find_global_maximum
