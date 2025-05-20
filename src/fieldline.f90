module fieldline_mod
    use constants, only: dp, pi
    implicit none

    type :: fieldline_t
        real(dp) :: theta_0
        real(dp) :: phi_0
        real(dp) :: iota
        real(dp) :: phi_max(2)
        real(dp) :: B_max(2)

        real(dp) :: eta_b
        real(dp) :: I_hat
    end type fieldline_t

contains

    subroutine make_flock_of_fieldlines(fieldlines, theta_0, iota)
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        real(dp), dimension(:), intent(in) :: theta_0
        real(dp), intent(in) :: iota

        fieldlines(:)%theta_0 = theta_0
        fieldlines(:)%iota = iota
    end subroutine make_flock_of_fieldlines

    subroutine set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, &
                                                   fieldlines)
        use field_base, only: field_t

        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta_mode, phi_mode
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines

        real(dp) :: chi_min_over_N

        call guess_chi_min_over_N(field, chi_min_over_N)
        fieldlines%phi_0 = theta_mode/phi_mode*fieldlines%theta_0 - chi_min_over_N
    end subroutine set_fieldline_phi_0_to_mode_minimum

    subroutine guess_chi_min_over_N(field, chi_min_over_N)
        use field_base, only: field_t
        use find_extrema, only: find_local_minima

        class(field_t), intent(in) :: field
        real(dp), intent(out) :: chi_min_over_N

        ! chi = M*theta - N*phi
        ! as f~f(chi) = sum c_j*cos(j*chi) with 1<j<j_max
        ! periodic at least after 2pi -> minimum must be in e.g [0, 3pi]
        real(dp), dimension(2), parameter :: interval = (/0.0_dp, 3.0_dp*pi/)
        real(dp) :: location(1)

        call find_local_minima(estimate_B_mod_of_chi_over_N, interval, location)
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
                                           n_steps)
        use find_extrema, only: find_local_maxima
        use field_base, only: field_t
        class(field_t), intent(in) :: field
        type(fieldline_t), intent(inout) :: fieldline
        real(dp), intent(in) :: interval(2)
        integer, intent(in), optional :: n_steps

        call find_local_maxima(B_mod_along_fieldline, interval, &
                               fieldline%phi_max, n_steps_in=n_steps)

        call B_mod_along_fieldline(fieldline%phi_max, fieldline%B_max)

    contains
        subroutine B_mod_along_fieldline(phi, B_mod)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(phi, 1)) :: theta
            integer :: idx

            theta = (phi - fieldline%phi_0)*fieldline%iota + fieldline%theta_0

            do idx = 1, size(phi, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_fieldline
    end subroutine find_maxima_along_fieldline

end module fieldline_mod
