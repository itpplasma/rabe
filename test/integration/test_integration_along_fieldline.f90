program test_fieldline
    use constants, only: dp, pi
    use utils, only: is_same
    use mock_field, only: mock_field_t

    implicit none

    real(dp), parameter :: theta_mode = 2.0_dp, phi_mode = -4.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, B_amplitude = -0.5_dp
    type(mock_field_t) :: field

    call field%mock_field_init(theta_mode, phi_mode, B_0, B_amplitude)
    call test_get_fieldlines()

contains

    subroutine test_get_fieldlines()
        use fieldline_mod, only: fieldline_t
        use fieldline_mod, only: set_fieldline_phi_0_to_mode_minimum
        use fieldline_mod, only: find_maxima_along_fieldline
        use utils, only: linspace

        real(dp), parameter :: reltol = 1e-2
        real(dp), parameter :: stor = 0.5_dp
        integer, parameter :: n_fieldlines = 10, n_maxima = 2

        real(dp), dimension(n_fieldlines) :: theta_0
        type(fieldline_t), dimension(n_fieldlines) :: fieldlines

        integer :: current
        real(dp) :: interval(2)
        real(dp) :: phi_max(n_maxima)

        call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, theta_0)
        fieldlines(:)%theta_0 = theta_0(:)
        fieldlines(:)%iota = -1.0_dp

        call set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, &
                                                 fieldlines)

        do current = 1, n_fieldlines
            allocate (fieldlines(current)%phi_max(n_maxima))
            interval = (/0.0_dp, 2*pi/) + fieldlines(current)%phi_0
            call find_maxima_along_fieldline(field, fieldlines(current), &
                                             interval, fieldlines(current)%phi_max)
            phi_max = (/0.5*pi, 1.5*pi/) + fieldlines(current)%phi_0
            if (is_same(phi_max, fieldlines(current)%phi_max, reltol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_get_fieldlines failed: phi_max"
                print *, "found: ", fieldlines(current)%phi_max
                print *, "expected: ", phi_max
                error stop
            end if
        end do

    end subroutine test_get_fieldlines

end program test_fieldline
