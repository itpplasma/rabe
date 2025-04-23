program test_fieldline
    use constants, only: dp, pi

    implicit none

    real(dp) :: retol = 1e-2

    call test_find_maxima_along_fieldline()

contains

    subroutine test_find_maxima_along_fieldline()
        use neo_field, only: neo_field_t
        use fieldline, only: find_maxima_along_fieldline

        character(len=*), parameter :: bc_filename = "input/single_mode_m_2_n_minus4.bc"
        real(dp), parameter :: theta_mode = 2.0_dp, phi_mode = -4.0_dp
        real(dp), parameter :: alpha_max = pi !cos(m*theta - n*nfp*phi) = cos(alpha)
        real(dp), dimension(2), parameter :: interval = (/0.0_dp, 2.0_dp*pi/)

        type(neo_field_t) :: field
        real(dp) :: stor(4), theta_0(4), iota(4)
        real(dp) :: found_phi(2), analytic_phi(2)
        integer :: idx

        stor = (/0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp/)
        theta_0 = (/3.0_dp/4.0_dp*pi, pi, -pi, -pi/)
        iota = (/1.00_dp, -3.00_dp, 1.00_dp, -3.00_dp/)

        call field%neo_field_init(bc_filename, stor(1))
        do idx = 1, 4
            call field%neo_change_stor(stor(idx))
            call find_maxima_along_fieldline(field, &
                                             iota(idx), &
                                             theta_0(idx), &
                                             interval, &
                                             found_phi)
            call find_analytic_maxima_along_fieldline(theta_mode, &
                                                      phi_mode, &
                                                      alpha_max, &
                                                      iota(idx), &
                                                      theta_0(idx), &
                                                      analytic_phi)
            if (any(abs(found_phi/analytic_phi - 1) > retol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_find_maxima_along_fieldline failed: phi at maxima"
                print *, "found: ", found_phi
                print *, "analytic: ", analytic_phi
                error stop
            end if
        end do
    end subroutine test_find_maxima_along_fieldline

    subroutine find_analytic_maxima_along_fieldline(theta_mode, &
                                                    phi_mode, &
                                                    alpha_max, &
                                                    iota, &
                                                    theta_0, &
                                                    phi)
        real(dp), intent(in) :: theta_mode, phi_mode, alpha_max, iota, theta_0
        real(dp), dimension(:), intent(out) :: phi

        integer :: current_maximum
        real(dp) :: offset, shift_value, first_valid_phi_max, shift

        offset = (alpha_max - theta_mode*theta_0)/(iota*theta_mode - phi_mode)
        shift_value = 2.0_dp*pi/(iota*theta_mode - phi_mode)

        if (offset < 0.0_dp) then
            shift = abs(int(offset/shift_value)) + 1.0_dp
            shift_value = abs(shift_value)
            first_valid_phi_max = offset + shift_value*shift
        else !(offset > 0.0_dp)
            shift = int(offset/shift_value)
            shift_value = abs(shift_value)
            first_valid_phi_max = offset - shift*shift_value
        end if

        shift = 0.0_dp
        do current_maximum = 1, size(phi, 1)
            phi(current_maximum) = first_valid_phi_max + shift_value*shift
            shift = shift + 1.0_dp
        end do

    end subroutine find_analytic_maxima_along_fieldline

end program test_fieldline
