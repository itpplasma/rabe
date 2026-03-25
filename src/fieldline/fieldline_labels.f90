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

    function suspect_omnigenous_origin_not_minimum(field, N_tor, M_pol, phi_tol)
        use find_extrema, only: find_local_minima
        use find_extrema, only: find_local_maxima
        use find_extrema, only: find_global_extrema
        use field_base, only: field_t
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

        class(field_t), intent(in) :: field
        real(dp), intent(in) :: N_tor, M_pol
        real(dp), intent(in), optional :: phi_tol
        logical :: suspect_omnigenous_origin_not_minimum

        real(dp), dimension(2) :: interval
        real(dp) :: tol
        integer, parameter :: n_min = 10
        real(dp) :: chis(n_min)
        integer :: n
        integer :: id_chi
        real(dp), parameter :: height_tol = 0.3_dp
        real(dp), dimension(2) :: extrema
        real(dp) :: B_min, B_max, B_range, B_at_origin, height
        logical :: origin_is_minimum, origin_is_maximum

        if (present(phi_tol)) then
            tol = phi_tol*abs(M_pol*pi - N_tor)
        else
            tol = 3.0_dp*1e-2
        end if

        suspect_omnigenous_origin_not_minimum = .false.
        origin_is_maximum = .false.
        origin_is_minimum = .false.

        interval = [-pi, pi]
        call find_local_minima(B_mod_along_pi_line, interval, chis, tol)
        n = count(.not. ieee_is_nan(chis))
        do id_chi = 1, n
            if (is_multiple_of_2pi(chis(id_chi), tol)) then
                origin_is_minimum = .true.
                exit
            end if
        end do

        if (.not. origin_is_minimum) then
            call find_local_maxima(B_mod_along_pi_line, interval, chis, tol)
            n = count(.not. ieee_is_nan(chis))
            do id_chi = 1, n
                if (is_multiple_of_2pi(chis(id_chi), tol)) then
                    origin_is_maximum = .true.
                    exit
                end if
            end do
        end if

        if (.not. origin_is_minimum .and. .not. origin_is_maximum) then
            print *, "warning: origin is neither local minimum nor maximum!"
            print *, "Stellarator symmetry is violated!"
            error stop
        end if

        extrema = find_global_extrema(B_mod_along_pi_line, interval, tol)
        B_min = extrema(1)
        B_max = extrema(2)
        B_range = B_max - B_min
        call field%compute_B_mod(0.0_dp, 0.0_dp, B_at_origin)
        height = B_at_origin - B_min
        if (height/B_range > height_tol) then
            print *, "Detected that origin in provided field is"
            print *, "a significant local maximum of height > "
            print *, height_tol*100.0_dp, "% of the total B range!"
            suspect_omnigenous_origin_not_minimum = .true.
        end if

    contains

        !> we go along the theta = pi*phi line so that
        !> chi = M_pol*theta - N_tor*phi can be inverted
        subroutine B_mod_along_pi_line(chi, B_mod)
            real(dp), dimension(:), intent(in) :: chi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(chi, 1)) :: phi, theta
            integer :: idx

            !> as M_pol and N_tor are whole numbers that must no both be zero
            !> M_pol*pi - N_tor should never zero
            if (abs(M_pol*pi - N_tor) < 1e-8) then
                print *, "Error: (M_pol*pi - N_tor) must not be (close) zero."
                print *, "abs(M_pol*iota - N_tor) = ", abs(M_pol*pi - N_tor)
                error stop
            end if

            phi = chi/(M_pol*pi - N_tor)
            theta = pi*phi

            do idx = 1, size(chi, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_pi_line

    end function suspect_omnigenous_origin_not_minimum

    function is_multiple_of_2pi(angle, tol_in)
        real(dp), intent(in) :: angle
        real(dp), intent(in), optional :: tol_in
        logical :: is_multiple_of_2pi

        real(dp) :: tol
        real(dp) :: remainder
        if (present(tol_in)) then
            tol = tol_in
        else
            tol = 3.0_dp*1e-2
        end if

        remainder = abs(mod(angle, 2.0_dp*pi))
        is_multiple_of_2pi = remainder < tol .or. abs(remainder - 2.0*pi) < tol
    end function is_multiple_of_2pi

end module fieldline_labels
