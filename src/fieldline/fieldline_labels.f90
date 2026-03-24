module fieldline_labels
    use constants, only: dp, pi

    implicit none

contains

    subroutine get_labels(max_n_fieldlines, iota, M_pol, N_tor, nfp, &
                          xi_0, approx_iota)
        use diophantine, only: rational_approx
        use utils, only: linspace
        integer, intent(in) :: max_n_fieldlines
        real(dp), intent(in) :: iota
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), dimension(:), allocatable, intent(out) :: xi_0
        real(dp), intent(out) :: approx_iota

        real(dp) :: iota_p, iota_p_approx
        integer :: p, q
        integer :: n_fieldlines

        iota_p = calc_iota_p(iota, M_pol, N_tor, nfp)
        !> The symmetry points xi_0=pi and xi_0=pi-iota_p have to be either
        !> part of the label grid or lie symmetric between two labels. Only then
        !> are the sampled points of a symmetric function themselvese symmetric
        !> in respect to those points. The grid is automatic symmetric in respect
        !> to pi if its equidistant between 0 and 2pi. Additionally, if one finds an
        !> approximated iota so that iota_p/2pi is rational p/q, and takes q as
        !> number of points excluding the endpoint 2pi, then
        !> $$
        !> dxi_0 = 2pi/q = iota_p/p
        !> $$
        !> and pi-iota_p is a whole number of steps away from pi and therefore
        !> also either part of the grid or symmetric between two labels.
        call rational_approx(iota_p/(2.0_dp*pi), max_n_fieldlines, p, q)
        iota_p_approx = 2.0_dp*pi*p/q
        approx_iota = calc_iota(iota_p_approx, M_pol, N_tor, nfp)

        !> Any multiple of q would work, but we want to use the largest one that
        !> is smaller than max_n_fieldlines to get the best resolution.
        n_fieldlines = q*(max_n_fieldlines/q)
        allocate (xi_0(n_fieldlines))

        call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, &
                      xi_0, include_endpoint=.false.)

    end subroutine get_labels

    function calc_iota_p(iota, M_pol, N_tor, nfp) result(iota_p)
        real(dp), intent(in) :: iota
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp) :: iota_p

        if (abs(iota*M_pol - N_tor) < 1.0e-10_dp) then
            error stop "iota*M_pol - N_tor must not be zero (resonant)!"
        end if
        if (M_pol**2.0_dp + N_tor**2.0_dp < 1.0e-10_dp) then
            error stop "M_pol**2 + N_tor**2 must not be zero!"
        end if
        iota_p = sign(pi, iota*M_pol - N_tor)/(M_pol**2.0_dp + N_tor**2.0_dp)* &
                 (M_pol + &
                  nfp*(N_tor*iota + M_pol)/(iota*M_pol - N_tor))

    end function calc_iota_p

    function calc_iota(iota_p, M_pol, N_tor, nfp) result(iota)
        real(dp), intent(in) :: iota_p
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp) :: iota

        real(dp) :: z, sigma

        sigma = 1.0_dp
        z = iota_p*(M_pol**2.0_dp + N_tor**2.0_dp)/pi*sigma - M_pol
        iota = (M_pol*nfp + N_tor*z)/(M_pol*z - N_tor*nfp)

        if (iota*M_pol - N_tor < 0.0_dp) then
            sigma = -1.0_dp
            z = iota_p*(M_pol**2.0_dp + N_tor**2.0_dp)/pi*sigma - M_pol
            iota = (M_pol*nfp + N_tor*z)/(M_pol*z - N_tor*nfp)
        end if

    end function calc_iota

    subroutine set_fieldline_labels_along_chi_min(field, M_pol, N_tor, nfp, &
                                                  fieldlines, phi_tol)
        use field_base, only: field_t
        use fieldline_mod, only: fieldline_t
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), intent(in), optional :: phi_tol
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines

        real(dp) :: chi_min, tol

        call guess_chi_min(field, chi_min, N_tor, M_pol, phi_tol)

        if (present(phi_tol)) then
            tol = phi_tol*3.0_dp
        else
            tol = 3.0_dp*1e-2
        end if
        if (not_multiple_of_2pi(chi_min, tol)) then
            print *, "error: found chi_min is not multiple of 2pi"
            print *, "chi_min: ", chi_min/pi, "[pi]"
            print *, "The minima contour of the ideal omnigenous configuration"
            print *, "must pass through (theta=0,phi=0)!"
            error stop
        end if

        fieldlines%theta_0 = N_tor*fieldlines%xi_0/nfp
        fieldlines%phi_0 = M_pol*fieldlines%xi_0/nfp
    end subroutine set_fieldline_labels_along_chi_min

    subroutine guess_chi_min(field, chi_min, N_tor, M_pol, tol)
        use find_extrema, only: find_local_minima
        use field_base, only: field_t

        class(field_t), intent(in) :: field
        real(dp), intent(out) :: chi_min
        real(dp), intent(in) :: N_tor, M_pol
        real(dp), intent(in), optional :: tol

        ! chi = M*theta - N*phi
        ! as f~f(chi) = sum c_j*cos(j*chi) with 1<j<j_max
        ! periodic at least after 2pi -> minimum must be in e.g [-pi, 2pi]
        real(dp), dimension(2), parameter :: interval = [0.0_dp, 3.0_dp*pi]
        real(dp) :: location(1)

        ! If f(theta,phi) approx f(chi = M*theta - N*phi) one can estiamte
        ! f(chi/N) by choosing 1 specific theta-phi combination for that
        ! chi value e.g.
        ! - phi=-chi/N and theta=0 or
        ! - phi=0 and theta=chi/M

        if (nint(N_tor) /= 0) then
            call find_local_minima(B_mod_along_phi_axis, interval, location, tol)
        elseif (nint(M_pol) /= 0) then
            call find_local_minima(B_mod_along_theta_axis, interval, location, tol)
        else
            print *, "error in guess_chi_min: M_pol=N_tor=0"
            print *, "M_pol and N_tor must not be both zero!"
            error stop
        end if

        chi_min = location(1)

    contains

        subroutine B_mod_along_phi_axis(chi, B_mod)
            real(dp), dimension(:), intent(in) :: chi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(chi, 1)) :: phi, theta
            integer :: idx

            phi = -chi/N_tor
            theta = 0.0_dp

            do idx = 1, size(chi, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_phi_axis

        subroutine B_mod_along_theta_axis(chi, B_mod)
            real(dp), dimension(:), intent(in) :: chi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(chi, 1)) :: phi, theta
            integer :: idx

            theta = chi/M_pol
            phi = 0.0_dp

            do idx = 1, size(chi, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_theta_axis

    end subroutine guess_chi_min

    function not_multiple_of_2pi(angle, tol)
        real(dp), intent(in) :: angle
        real(dp), intent(in) :: tol
        logical :: not_multiple_of_2pi

        real(dp) :: remainder

        remainder = abs(mod(angle, 2.0_dp*pi))
        not_multiple_of_2pi = remainder > tol .and. abs(remainder - 2.0*pi) > tol
    end function not_multiple_of_2pi

end module fieldline_labels
