module fieldline
    use constants, only: dp, pi
    implicit none
contains
    subroutine guess_alpha_at_minimum(field, alpha_at_min, M)
        use field_base, only: field_t
        use find_extrema, only: find_local_minima

        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M
        real(dp), intent(out) :: alpha_at_min

        real(dp), dimension(2), parameter :: interval = (/0.0_dp, 4.0_dp*pi/)
        real(dp) :: location(1)

        call find_local_minima(estimate_B_mod_of_alpha, interval, location)
        alpha_at_min = mod(location(1), 2*pi)

    contains

        subroutine estimate_B_mod_of_alpha(alpha, B_mod)
            real(dp), dimension(:), intent(in) :: alpha
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(alpha, 1)) :: phi, theta
            integer :: idx

            ! If f(theta,phi) approx f(alpha = M*theta - N*phi) one can estiamte
            ! f(alpha) by choosing 1 specific theta-phi combination for that
            ! alpha value e.g. phi=0 and theta=alpha/M
            phi = 0.0_dp
            theta = alpha/M

            do idx = 1, size(alpha, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine estimate_B_mod_of_alpha

    end subroutine guess_alpha_at_minimum

    subroutine find_maxima_along_fieldline(field, iota, theta_0, phi_at_max)
        use find_extrema, only: find_local_maxima
        use field_base, only: field_t

        real(dp), parameter :: M = 1.0_dp
        real(dp), parameter :: nfp = 4.0_dp, scan_n_periods = 2.0_dp
        real(dp), dimension(2), parameter :: interval = (/0.0_dp, &
                                                          2.0_dp*pi/nfp*scan_n_periods/)

        class(field_t), intent(in) :: field
        real(dp), intent(in) :: iota, theta_0
        real(dp), dimension(:) :: phi_at_max

        call find_local_maxima(B_mod_along_fieldline, interval, phi_at_max)

    contains
        subroutine B_mod_along_fieldline(phi, B_mod)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(phi, 1)) :: theta
            integer :: idx

            theta = phi*iota + theta_0

            do idx = 1, size(phi, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_fieldline
    end subroutine find_maxima_along_fieldline
end module fieldline
