module make_fieldline
    use constants, only: dp, pi, eps
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t

    implicit none

    type :: surface_average_t
        real(dp) :: normalization
        real(dp) :: B_squared
        real(dp) :: lambda_b
    end type surface_average_t

contains

    subroutine make_flock_of_fieldlines(fieldlines, theta_0, iota, &
                                        field, M_pol, N_tor, phi_tol)
        use fieldline_integrals, only: calc_fieldline_integrals
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        real(dp), dimension(:), intent(in) :: theta_0
        real(dp), intent(in) :: iota
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor
        real(dp), intent(in), optional :: phi_tol

        real(dp) :: interval(2)
        real(dp) :: average_delta_aspect_ratio
        integer :: n_fieldlines
        integer :: current

        n_fieldlines = size(fieldlines)

        fieldlines%theta_0 = theta_0
        fieldlines%iota = iota

        call set_fieldline_phi_0_to_mode_minimum(field, M_pol, N_tor, fieldlines, &
                                                 phi_tol)

        fieldlines%iota_p = iota*(fieldlines(1)%theta_0*M_pol &
                                  - fieldlines(2)%phi_0*N_tor)/(N_tor - iota*M_pol)

        do current = 1, n_fieldlines
            interval = (/-1.5_dp*pi, 1.5_dp*pi/) + fieldlines(current)%phi_0
            call find_maxima_along_fieldline(field, fieldlines(current), &
                                             interval, phi_tol)
        end do

        fieldlines%eta_b = (1.0_dp - eps)/get_global_B_max(fieldlines)
        fieldlines%delta_eta = 1.0_dp/fieldlines(:)%B_max(1) - fieldlines(:)%eta_b

        do current = 1, n_fieldlines
            call calc_fieldline_integrals(field, fieldlines(current))
        end do

        ! I_ref can be chosen to be any e.g. I = I_1
        ! (I_ref/I_j)**0.5 - 1 = (max(I_j)/I_j)**0.5 -1 =
        ! ((I+delta)/(I+delta_j))**0.5 -1 ~ 0.5*(delta/I - delta_j/I)
        ! and the result in linear order only differs by a constant delta/I
        ! which does not enter the offset formula
        fieldlines(:)%delta_aspect_ratio = sqrt( &
                                fieldlines(1)%integral_lambda_b_over_B_squared/ &
                                fieldlines(:)%integral_lambda_b_over_B_squared &
                                           ) - 1
        ! average of delta_aspect ratio also does not enter offset formula
        ! can be set it to zero
        average_delta_aspect_ratio = sum(fieldlines(:)%delta_aspect_ratio)/n_fieldlines
        fieldlines(:)%delta_aspect_ratio = fieldlines(:)%delta_aspect_ratio - &
                                           average_delta_aspect_ratio

    end subroutine make_flock_of_fieldlines

    subroutine set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, &
                                                   fieldlines, phi_tol)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta_mode, phi_mode
        real(dp), intent(in), optional :: phi_tol
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines

        real(dp) :: chi_min_over_N

        call guess_chi_min_over_N(field, chi_min_over_N, phi_tol)
        fieldlines%phi_0 = theta_mode/phi_mode*fieldlines%theta_0 - chi_min_over_N
    end subroutine set_fieldline_phi_0_to_mode_minimum

    subroutine guess_chi_min_over_N(field, chi_min_over_N, phi_tol)
        use find_extrema, only: find_local_minima

        class(field_t), intent(in) :: field
        real(dp), intent(in), optional :: phi_tol
        real(dp), intent(out) :: chi_min_over_N

        ! chi = M*theta - N*phi
        ! as f~f(chi) = sum c_j*cos(j*chi) with 1<j<j_max
        ! periodic at least after 2pi -> minimum must be in e.g [0, 3pi]
        real(dp), dimension(2), parameter :: interval = (/0.0_dp, 3.0_dp*pi/)
        real(dp) :: location(1)

        call find_local_minima(estimate_B_mod_of_chi_over_N, interval, location, &
                               phi_tol)
        chi_min_over_N = location(1)

    contains

        subroutine estimate_B_mod_of_chi_over_N(chi_over_N, B_mod)
            real(dp), dimension(:), intent(in) :: chi_over_N
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(chi_over_N, 1)) :: phi, theta
            integer :: idx

            ! If f(theta,phi) approx f(chi = M*theta - N*phi) one can estiamte
            ! f(chi/N) by choosing 1 specific theta-phi combination for that
            ! chi value e.g. phi=-chi/N and theta=0
            phi = -chi_over_N
            theta = 0.0_dp

            do idx = 1, size(chi_over_N, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine estimate_B_mod_of_chi_over_N

    end subroutine guess_chi_min_over_N

    subroutine find_maxima_along_fieldline(field, &
                                           fieldline, &
                                           interval, &
                                           phi_tol)
        use find_extrema, only: find_local_maxima
        class(field_t), intent(in) :: field
        type(fieldline_t), intent(inout) :: fieldline
        real(dp), intent(in) :: interval(2)
        real(dp), intent(in), optional :: phi_tol

        call find_local_maxima(B_mod_along_fieldline, interval, &
                               fieldline%phi_max, phi_tol)

        call B_mod_along_fieldline(fieldline%phi_max, fieldline%B_max)

    contains
        subroutine B_mod_along_fieldline(phi, B_mod)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(phi)) :: theta
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_fieldline
    end subroutine find_maxima_along_fieldline

    function get_global_B_max(fieldlines) result(global_B_max)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: global_B_max

        integer :: current
        real(dp) :: B_max_1, B_max_2

        B_max_1 = maxval(fieldlines(:)%B_max(1))
        B_max_2 = maxval(fieldlines(:)%B_max(2))
        global_B_max = max(B_max_1, B_max_2)
    end function get_global_B_max

    subroutine calc_surface_averages(fieldlines, surface_average)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(surface_average_t), intent(out) :: surface_average

        real(dp), dimension(size(fieldlines)) :: well_lengths

        well_lengths = fieldlines%phi_max(2) - fieldlines%phi_max(1)
        surface_average%normalization = sum(fieldlines%integral_one_over_B_squared)
        surface_average%B_squared = sum(well_lengths)/surface_average%normalization
        surface_average%lambda_b = sum(fieldlines%integral_lambda_b_over_B_squared)/ &
                                   surface_average%normalization
    end subroutine calc_surface_averages

end module make_fieldline
