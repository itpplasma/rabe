module fieldline_mod
    use constants, only: dp, pi, eps
    use field_base, only: field_t

    implicit none

    type :: fieldline_t
        real(dp) :: theta_0
        real(dp) :: phi_0
        real(dp) :: iota
        real(dp) :: phi_max(2)
        real(dp) :: B_max(2)
        real(dp) :: iota_p

        real(dp) :: eta_b
        real(dp) :: delta_eta
        real(dp) :: well_average_lambda
        real(dp) :: delta_aspect_ratio
        real(dp) :: well_average_lambda_b
        real(dp) :: radial_drift
    contains
        generic :: get_theta => get_theta_scalar, get_theta_array
        procedure, private :: get_theta_scalar
        procedure, private :: get_theta_array
    end type fieldline_t

contains

    function get_theta_scalar(self, phi) result(theta)
        class(fieldline_t), intent(in) :: self
        real(dp) :: phi

        real(dp) :: theta

        theta = (phi - self%phi_0)*self%iota + self%theta_0
    end function get_theta_scalar

    function get_theta_array(self, phi) result(theta)
        class(fieldline_t), intent(in) :: self
        real(dp), dimension(:) :: phi

        real(dp), dimension(size(phi)) :: theta

        theta = (phi - self%phi_0)*self%iota + self%theta_0
    end function get_theta_array

    subroutine make_flock_of_fieldlines(fieldlines, theta_0, iota, &
                                        field, M_pol, N_tor, phi_tol)
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        real(dp), dimension(:), intent(in) :: theta_0
        real(dp), intent(in) :: iota
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor
        real(dp), intent(in), optional :: phi_tol

        real(dp) :: interval(2)
        integer :: current

        fieldlines%theta_0 = theta_0
        fieldlines%iota = iota

        call set_fieldline_phi_0_to_mode_minimum(field, M_pol, N_tor, fieldlines, &
                                                 phi_tol)

        fieldlines%iota_p = iota*(fieldlines(1)%theta_0*M_pol &
                                  - fieldlines(2)%phi_0*N_tor)/(N_tor - iota*M_pol)

        do current = 1, size(fieldlines)
            interval = (/-1.5_dp*pi, 1.5_dp*pi/) + fieldlines(current)%phi_0
            call find_maxima_along_fieldline(field, fieldlines(current), &
                                             interval, phi_tol)
        end do

        fieldlines%eta_b = (1.0_dp - eps)/get_global_B_max(fieldlines)
        fieldlines%delta_eta = 1.0_dp/fieldlines(:)%B_max(1) - fieldlines(:)%eta_b

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

end module fieldline_mod
