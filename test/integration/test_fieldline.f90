program test_fieldline
    use constants, only: dp, pi
    use neo_field, only: neo_field_t
    use utils, only: is_same

    implicit none

    character(len=*), parameter :: bc_filename = "input/single_mode_m_2_n_minus4.bc"
    real(dp), parameter :: theta_mode = 2.0_dp, phi_mode = -4.0_dp
    !The minimum/maximum alpha of a single mode field
    !-cos(alpha) with alpha = M*theta - N*phi
    real(dp), parameter :: alpha_max = pi
    real(dp), parameter :: alpha_min = 0.0_dp

    type(neo_field_t) :: field

    call field%neo_field_init(bc_filename, 0.0_dp)
    call test_guess_alpha_at_minimum()
    call test_find_maxima_along_fieldline()
    call test_set_fieldline_labels_to_mode_minimum()
    call test_get_fieldlines()

contains

    subroutine test_guess_alpha_at_minimum()
        use fieldline_mod, only: guess_alpha_over_M_at_minimum

        real(dp), parameter :: retol = 1e-3

        real(dp) :: stor(4)
        real(dp) :: found_alpha_over_M_min, found_alpha_min
        integer :: idx

        stor = (/0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp/)
        do idx = 1, 4
            call field%neo_change_stor(stor(idx))
            call guess_alpha_over_M_at_minimum(field, found_alpha_over_M_min)
            found_alpha_min = mod(found_alpha_over_M_min*theta_mode + pi, 2*pi) - pi
            if (is_same(alpha_min, found_alpha_min, retol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_guess_alpah_at_minimum failed: alpha at minima"
                print *, "found: ", found_alpha_min
                print *, "analytic: ", alpha_min
                error stop
            end if
        end do
    end subroutine test_guess_alpha_at_minimum

    subroutine test_find_maxima_along_fieldline()
        use fieldline_mod, only: find_maxima_along_fieldline, fieldline_t

        real(dp), parameter :: retol = 1e-2, abstol = 0.0_dp
        integer, parameter :: n_maxima = 2, n_steps = 1000

        real(dp) :: stor(4), theta_0(4), phi_0(4), iota(4)
        real(dp) :: interval(2)
        type(fieldline_t) :: fieldline
        real(dp) :: found_phi(2), analytic_phi(2)
        integer :: idx

        stor = (/0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp/)
        theta_0 = (/3.0_dp/4.0_dp*pi, pi, -pi, -pi/)
        phi_0 = (/0.25_dp*pi, 1.0_dp/3.0_dp*pi, 1.25_dp*pi, 1.50_dp*pi/)
        iota = (/1.00_dp, -3.00_dp, 1.00_dp, -3.00_dp/)

        allocate (fieldline%phi_max(n_maxima))
        do idx = 1, 4
            call field%neo_change_stor(stor(idx))
            fieldline%theta_0 = theta_0(idx)
            fieldline%phi_0 = phi_0(idx)
            fieldline%iota = iota(idx)
            interval = (/0.0_dp, 2.0_dp*pi/) + phi_0(idx)
            call find_maxima_along_fieldline(field, &
                                             fieldline, &
                                             interval, &
                                             n_steps)
            call find_analytic_maxima_along_fieldline(theta_mode, &
                                                      phi_mode, &
                                                      alpha_max, &
                                                      fieldline, &
                                                      analytic_phi)
            if (is_same(analytic_phi, fieldline%phi_max, retol, abstol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_find_maxima_along_fieldline failed: phi at maxima"
                print *, "found: ", fieldline%phi_max
                print *, "analytic: ", analytic_phi
                error stop
            end if
        end do
    end subroutine test_find_maxima_along_fieldline

    subroutine find_analytic_maxima_along_fieldline(m, n, alpha_max, &
                                                    fieldline, phi)
        use fieldline_mod, only: fieldline_t
        real(dp), intent(in) :: m, n, alpha_max
        type(fieldline_t) :: fieldline
        real(dp), dimension(:), intent(out) :: phi

        integer :: current_maximum
        real(dp) :: offset, shift_value, first_valid_phi_max, shift

        offset = (alpha_max - m*fieldline%theta_0 + n*fieldline%phi_0)/ &
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
        use fieldline_mod, only: fieldline_t, set_fieldline_phi_0_to_mode_minimum
        use utils, only: linspace

        real(dp), parameter :: retol = 1e-8, abstol = 0.0_dp
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
            if (is_same(B_mod_min, B_mod, retol, abstol)) then
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

    subroutine test_get_fieldlines()
        use fieldline_mod, only: fieldline_t
        use fieldline_mod, only: set_fieldline_phi_0_to_mode_minimum
        use fieldline_mod, only: find_maxima_along_fieldline
        use utils, only: linspace

        real(dp), parameter :: reltol = 1e-2
        integer, parameter :: n_steps = 1000
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

        call field%neo_change_stor(stor)
        call set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, &
                                                 fieldlines)

        do current = 1, n_fieldlines
            allocate (fieldlines(current)%phi_max(n_maxima))
            interval = (/0.0_dp, 2*pi/) + fieldlines(current)%phi_0
            call find_maxima_along_fieldline(field, fieldlines(current), &
                                             interval, n_steps)
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
