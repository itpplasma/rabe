module fieldline_mod
    use constants, only: dp, pi
    implicit none

    type :: fieldline_t
        real(dp) :: theta_0
        real(dp) :: phi_0
        real(dp) :: iota
        real(dp), dimension(:), allocatable :: phi_max
    end type fieldline_t

contains

    subroutine set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, &
                                                   fieldlines)
        use field_base, only: field_t

        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta_mode, phi_mode
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines

        real(dp) :: alpha_over_M_min

        call guess_alpha_over_M_at_minimum(field, alpha_over_M_min)
        fieldlines%phi_0 = theta_mode/phi_mode*(fieldlines%theta_0 - alpha_over_M_min)
    end subroutine set_fieldline_phi_0_to_mode_minimum

    subroutine guess_alpha_over_M_at_minimum(field, alpha_over_M_min)
        use field_base, only: field_t
        use find_extrema, only: find_local_minima

        class(field_t), intent(in) :: field
        real(dp), intent(out) :: alpha_over_M_min

        ! alpha = M*theta - N*phi
        ! as f~f(alpha) = sum c_j*cos(j*alpha) with 1<j<j_max
        ! periodic at least after 2pi -> minimum must be in e.g [0, 3pi]
        real(dp), dimension(2), parameter :: interval = (/0.0_dp, 3.0_dp*pi/)
        real(dp) :: location(1)

        call find_local_minima(estimate_B_mod_of_alpha_over_M, interval, location)
        alpha_over_M_min = location(1)

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

    subroutine find_maxima_along_fieldline(field, fieldline, interval, phi_at_max)
        use find_extrema, only: find_local_maxima
        use field_base, only: field_t
        class(field_t), intent(in) :: field
        type(fieldline_t), intent(in) :: fieldline
        real(dp), intent(in) :: interval(2)
        real(dp), dimension(:) :: phi_at_max

        call find_local_maxima(B_mod_along_fieldline, interval, phi_at_max)

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
