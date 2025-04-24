module fieldline
    use constants, only: dp, pi
    implicit none

contains
    subroutine guess_alpha_over_M_at_minimum(field, alpha_over_M_min)
        use field_base, only: field_t
        use find_extrema, only: find_local_minima

        class(field_t), intent(in) :: field
        real(dp), intent(out) :: alpha_over_M_min

        real(dp), dimension(2), parameter :: interval = (/-pi/4.0_dp, 1.75_dp*pi/)
        real(dp) :: location(1)

        call find_local_minima(estimate_B_mod_of_alpha_over_M, interval, location)
        alpha_over_M_min = mod(location(1), 2*pi)

    contains

        subroutine estimate_B_mod_of_alpha_over_M(alpha_over_M, B_mod)
            real(dp), dimension(:), intent(in) :: alpha_over_M
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(alpha_over_M, 1)) :: phi, theta
            integer :: idx

            ! If f(theta,phi) approx f(alpha = M*theta - N*phi) one can estiamte
            ! f(alpha/M) by choosing 1 specific theta-phi combination for that
            ! alpha value e.g. phi=0 and theta=alpha/M
            phi = 0.0_dp
            theta = alpha_over_M

            do idx = 1, size(alpha_over_M, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine estimate_B_mod_of_alpha_over_M

    end subroutine guess_alpha_over_M_at_minimum

    subroutine find_maxima_along_fieldline(field, iota, theta_0, interval, phi_at_max)
        use find_extrema, only: find_local_maxima
        use field_base, only: field_t
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: iota, theta_0, interval(2)
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
