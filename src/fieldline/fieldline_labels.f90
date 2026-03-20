module fieldline_labels
    use constants, only: dp, pi

    implicit none

contains

    subroutine get_theta_0(max_n_fieldlines, iota, M_pol, N_tor, nfp, &
                           theta_0, approx_iota)
        use diophantine, only: rational_approx
        use utils, only: linspace
        integer, intent(in) :: max_n_fieldlines
        real(dp), intent(in) :: iota
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), dimension(:), allocatable, intent(out) :: theta_0
        real(dp), intent(out) :: approx_iota

        real(dp) :: iota_p, iota_p_approx
        integer :: p, q
        integer :: n_fieldlines

        iota_p = calc_iota_p(iota, M_pol, N_tor, nfp)
        !> The symmetry points theta_0=pi and theta_0=pi-iota_p have to be either
        !> part of the label grid or lie symmetric between two labels. Only then
        !> are the sampled points of a symmetric function themselvese symmetric
        !> in respect to those points. The grid is automatic symmetric in respect
        !> to pi if its equidistant between 0 and 2pi. Additionally, if one finds an
        !> approximated iota so that iota_p/2pi is rational p/q, and takes q as
        !> number of points excluding the endpoint 2pi, then
        !> $$
        !> dtheta_0 = 2pi/q = iota_p/p
        !> $$
        !> and pi-iota_p is a whole number of steps away from pi and therefore
        !> also either part of the grid or symmetric between two labels.
        call rational_approx(iota_p/(2.0_dp*pi), max_n_fieldlines, p, q)
        iota_p_approx = 2.0_dp*pi*p/q
        approx_iota = calc_iota(iota_p_approx, M_pol, N_tor, nfp)

        !> Any multiple of q would work, but we want to use the largest one that
        !> is smaller than max_n_fieldlines to get the best resolution.
        n_fieldlines = q*(max_n_fieldlines/q)
        allocate (theta_0(n_fieldlines))

        call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, &
                      theta_0, include_endpoint=.false.)

    end subroutine get_theta_0

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

end module fieldline_labels
