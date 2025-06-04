program test_fieldline
    use constants, only: dp, pi
    use neo_field, only: neo_field_t
    use utils, only: not_same

    implicit none

    character(len=*), parameter :: bc_filename = "input/single_mode_m_2_n_minus4.bc"
    real(dp), parameter :: theta_mode = 2.0_dp, phi_mode = -4.0_dp
    !The minimum/maximum chi of a single mode field
    !-cos(chi) with chi = M*theta - N*phi
    real(dp), parameter :: chi_max = pi
    real(dp), parameter :: chi_min = 0.0_dp

    type(neo_field_t) :: field

    call field%neo_field_init(bc_filename, 0.0_dp)
    call test_guess_chi_at_minimum()
    call test_find_maxima_along_fieldline()
    call test_set_fieldline_labels_to_mode_minimum()

contains

    subroutine test_guess_chi_at_minimum()
        use make_fieldline, only: guess_chi_min_over_N

        real(dp), parameter :: retol = 1e-2*0.5_dp*abs(phi_mode)

        real(dp) :: stor(4)
        real(dp) :: found_chi_min_over_N, found_chi_min
        integer :: idx

        stor = (/0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp/)
        do idx = 1, 4
            call field%neo_change_stor(stor(idx))
            call guess_chi_min_over_N(field, found_chi_min_over_N)
            found_chi_min = mod(found_chi_min_over_N*phi_mode - pi, 2*pi) + pi
            if (not_same(chi_min, found_chi_min, retol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_guess_chi_min failed: chi at minima"
                print *, "found: ", found_chi_min
                print *, "analytic: ", chi_min
                error stop
            end if
        end do
    end subroutine test_guess_chi_at_minimum

    subroutine test_find_maxima_along_fieldline()
        use make_fieldline, only: find_maxima_along_fieldline
        use fieldline_mod, only: fieldline_t

        real(dp), parameter :: phi_tol = 1e-3
        real(dp), parameter :: retol = 0.0_dp, abstol = phi_tol

        real(dp) :: stor(4), theta_0(4), phi_0(4), iota(4)
        real(dp) :: interval(2)
        type(fieldline_t) :: fieldline
        real(dp) :: found_phi(2), analytic_phi(2)
        integer :: idx

        stor = (/0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp/)
        theta_0 = (/3.0_dp/4.0_dp*pi, pi, -pi, -pi/)
        phi_0 = (/0.25_dp*pi, 1.0_dp/3.0_dp*pi, 1.25_dp*pi, 1.50_dp*pi/)
        iota = (/1.00_dp, -3.00_dp, 1.00_dp, -3.00_dp/)

        do idx = 1, 4
            call field%neo_change_stor(stor(idx))
            fieldline%theta_0 = theta_0(idx)
            fieldline%phi_0 = phi_0(idx)
            fieldline%iota = iota(idx)
            interval = (/0.0_dp, 2.0_dp*pi/) + phi_0(idx)
            call find_maxima_along_fieldline(field, &
                                             fieldline, &
                                             interval, &
                                             phi_tol)
            call find_analytic_maxima_along_fieldline(theta_mode, &
                                                      phi_mode, &
                                                      chi_max, &
                                                      fieldline, &
                                                      analytic_phi)
            if (not_same(analytic_phi, fieldline%phi_max, retol, abstol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_find_maxima_along_fieldline failed: phi at maxima"
                print *, "found: ", fieldline%phi_max
                print *, "analytic: ", analytic_phi
                error stop
            end if
        end do
    end subroutine test_find_maxima_along_fieldline

    subroutine find_analytic_maxima_along_fieldline(m, n, chi_max, &
                                                    fieldline, phi)
        use fieldline_mod, only: fieldline_t
        real(dp), intent(in) :: m, n, chi_max
        type(fieldline_t) :: fieldline
        real(dp), dimension(:), intent(out) :: phi

        integer :: current_maximum
        real(dp) :: offset, shift_value, first_valid_phi_max, shift

        offset = (chi_max - m*fieldline%theta_0 + n*fieldline%phi_0)/ &
                 (fieldline%iota*m - n)
        shift_value = 2.0_dp*pi/(fieldline%iota*m - n)

        if (offset < 0.0_dp) then
            shift = abs(int(offset/shift_value)) + 1.0_dp
            shift_value = abs(shift_value)
            first_valid_phi_max = offset + shift_value*shift + fieldline%phi_0
        else
            shift = abs(int(offset/shift_value))
            shift_value = abs(shift_value)
            first_valid_phi_max = offset - shift*shift_value + fieldline%phi_0
        end if

        shift = 0.0_dp
        do current_maximum = 1, size(phi, 1)
            phi(current_maximum) = first_valid_phi_max + shift_value*shift
            shift = shift + 1.0_dp
        end do

    end subroutine find_analytic_maxima_along_fieldline

    subroutine test_set_fieldline_labels_to_mode_minimum()
        use fieldline_mod, only: fieldline_t
        use make_fieldline, only: set_fieldline_phi_0_to_mode_minimum
        use utils, only: linspace

        real(dp), parameter :: retol = (1e-2*phi_mode)**2, abstol = 0.0_dp
        real(dp), parameter :: stor = 0.5_dp
        integer, parameter :: n_fieldlines = 10

        real(dp), dimension(n_fieldlines) :: theta_0
        type(fieldline_t), dimension(n_fieldlines) :: fieldlines

        real(dp) :: B_mod
        real(dp) :: B_mod_min = 0.7_dp
        integer :: current

        call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, theta_0)
        fieldlines(:)%theta_0 = theta_0(:)

        call field%neo_change_stor(stor)
        call set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, &
                                                 fieldlines)
        do current = 1, size(fieldlines)
            call field%compute_B_mod(fieldlines(current)%theta_0, &
                                     fieldlines(current)%phi_0, &
                                     B_mod)
            if (not_same(B_mod_min, B_mod, retol, abstol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_set_fieldline_labels_to_mode_minimum failed: B_mod"
                print *, "at theta_0", fieldlines(current)%theta_0
                print *, "at phi_0", fieldlines(current)%phi_0
                print *, "found: ", B_mod
                print *, "minimum: ", B_mod_min
                error stop
            end if
        end do
    end subroutine test_set_fieldline_labels_to_mode_minimum

end program test_fieldline
