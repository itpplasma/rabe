program test_fieldline
    use constants, only: dp, pi
    use neo_field, only: neo_field_t
    use utils, only: not_same

    implicit none

    character(len=*), parameter :: bc_filename = "input/single_mode_m_2_n_minus4.bc"
    real(dp), parameter :: M_pol = -2.0_dp, N_tor = 4.0_dp
    real(dp), parameter :: nfp = max(1.0_dp, abs(N_tor))
    !The minimum/maximum chi of a single mode field
    !-cos(chi) with chi = M*theta - N*phi
    real(dp), parameter :: chi_max = pi
    real(dp), parameter :: chi_min = 0.0_dp

    type(neo_field_t) :: field

    call field%neo_field_init(bc_filename, 0.0_dp)
    call test_find_maxima_along_fieldline()
    call test_set_fieldline_labels_to_mode_minimum()

contains

    subroutine test_find_maxima_along_fieldline()
        use make_fieldline, only: find_maxima_along_fieldline
        use make_fieldline, only: maxima_t
        use fieldline_mod, only: fieldline_t

        real(dp), parameter :: retol = 0.0_dp, abstol = 1e-6

        real(dp) :: stor(4), theta_0(4), phi_0(4), iota(4)
        real(dp) :: interval(2)
        type(fieldline_t) :: fieldline
        type(maxima_t) :: maxima
        real(dp) :: found_phi(2), analytic_phi(2)
        integer :: idx

        stor = [0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp]
        theta_0 = [3.0_dp/4.0_dp*pi, pi, -pi, -pi]
        phi_0 = [0.25_dp*pi, 1.0_dp/3.0_dp*pi, 1.25_dp*pi, 1.50_dp*pi]
        iota = [1.00_dp, -3.00_dp, 1.00_dp, -3.00_dp]

        do idx = 1, 4
            call field%neo_change_stor(stor(idx))
            fieldline%theta_0 = theta_0(idx)
            fieldline%phi_0 = phi_0(idx)
            fieldline%iota = iota(idx)
            interval = [0.0_dp, 4.1_dp*pi]/abs(M_pol*iota(idx) - N_tor) + phi_0(idx)
            call find_maxima_along_fieldline(field, &
                                             fieldline, &
                                             interval, &
                                             maxima)
            if (maxima%n == 2) then
                fieldline%phi_max = maxima%phi(1:2)
            else
                print *, "-------------------------------------------------------------"
                print *, "test_integration_along_fieldline failed: maxima"
                print *, "Found ", maxima%n, " maxima, expected 2!"
                print *, "maxima: ", maxima%phi(1:max(1, maxima%n))
                error stop
            end if

            call find_analytic_maxima_along_fieldline(M_pol, &
                                                      N_tor, &
                                                      chi_max, &
                                                      fieldline, &
                                                      analytic_phi)
            if (not_same(analytic_phi, fieldline%phi_max, retol, abstol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_find_maxima_along_fieldline failed: phi at maxima"
                print *, "case: ", idx
                print *, "found: ", fieldline%phi_max
                print *, "analytic: ", analytic_phi
                print *, "error: ", abs(fieldline%phi_max - analytic_phi)
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
        use fieldline_labels, only: suspect_omnigenous_origin_not_minimum
        use utils, only: linspace

        real(dp), parameter :: retol = (1e-2*N_tor)**2, abstol = 0.0_dp
        real(dp), parameter :: stor = 0.5_dp
        integer, parameter :: n_fieldlines = 10

        real(dp), dimension(n_fieldlines) :: xi_0
        type(fieldline_t), dimension(n_fieldlines) :: fieldlines

        real(dp) :: B_mod
        real(dp) :: B_mod_min = 0.7_dp
        integer :: current
        logical :: test_failed

        test_failed = .false.
        call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, xi_0)
        fieldlines(:)%xi_0 = xi_0(:)

        call field%neo_change_stor(stor)
        if (suspect_omnigenous_origin_not_minimum(field, M_pol, N_tor)) then
            print *, "error: The origin of the IDEAL omnigenous configuration"
            print *, "(theta=phi=0) must be a global and local minimum!"
            print *, "Origin of provided field suggests that this is not the case!"
            error stop
        end if
        fieldlines%theta_0 = N_tor*fieldlines%xi_0/nfp
        fieldlines%phi_0 = M_pol*fieldlines%xi_0/nfp
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
                test_failed = .true.
            end if
        end do

        if (test_failed) then
            error stop
        end if
    end subroutine test_set_fieldline_labels_to_mode_minimum

end program test_fieldline
