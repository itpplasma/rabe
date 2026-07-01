module boozmn_reader
    !> Load Boozer splines from a booz_xform boozmn NetCDF.
    !>
    !> Reads the Fourier harmonics bmnc_b, iota_b, buco_b, bvco_b, phi_b and the
    !> mode arrays ixm_b, ixn_b, evaluates Bmod on a uniform grid by Fourier
    !> summation, and hands the result to libneo's build_boozer_from_chartmap so
    !> the boozmn and chartmap paths share one spline backend.

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: load_boozer_from_boozmn

    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: GAUSS_CM2_PER_TM2 = 1.0e8_dp
    real(dp), parameter :: GAUSS_CM_PER_TM = 1.0e6_dp
    real(dp), parameter :: GAUSS_PER_T = 1.0e4_dp

contains

    subroutine load_boozer_from_boozmn(boozmn_file, nrho_in, ntheta_in, nzeta_in)
        use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
        use boozer_chartmap_types, only: boozer_chartmap_data_t
        use boozer_sub, only: build_boozer_from_chartmap

        character(len=*), intent(in) :: boozmn_file
        integer, intent(in), optional :: nrho_in, ntheta_in, nzeta_in

        type(boozer_chartmap_data_t) :: d
        integer :: nrho, ntheta, nzeta
        integer :: ncid, ns, nmn, nsurf, nfp_int, lasym
        integer :: ir, it, iz, mn, mn00, k
        integer, allocatable :: jlist(:), ixm(:), ixn(:)
        real(dp), allocatable :: iota_full(:), buco_full(:), bvco_full(:), phi_full(:)
        real(dp), allocatable :: bmnc_h(:, :), rmnc_h(:, :)
        real(dp), allocatable :: rho_half(:), s_half(:)
        real(dp), allocatable :: theta(:), zeta(:), bmnc_out(:, :)
        real(dp), allocatable :: bmod_geom(:, :, :)
        real(dp) :: torflux_si, angle, rmajor_m
        real(dp) :: rho_min, rho_max, s_min, s_max

        nrho = 30
        ntheta = 48
        nzeta = 96
        if (present(nrho_in)) nrho = nrho_in
        if (present(ntheta_in)) ntheta = ntheta_in
        if (present(nzeta_in)) nzeta = nzeta_in

        call nc_open(trim(boozmn_file), ncid)
        call nc_get(ncid, 'ns_b', ns)
        call nc_inq_dim(ncid, 'ixm_b', nmn)
        call nc_inq_dim(ncid, 'jlist', nsurf)

        allocate (jlist(nsurf), ixm(nmn), ixn(nmn))
        allocate (iota_full(ns), buco_full(ns), bvco_full(ns), phi_full(ns))
        allocate (bmnc_h(nmn, nsurf), rmnc_h(nmn, nsurf))

        call nc_get(ncid, 'nfp_b', nfp_int)
        call nc_get(ncid, 'lasym__logical__', lasym)
        if (lasym /= 0) then
            call nc_close(ncid)
            error stop "asymmetric (lasym=1) boozmn not yet supported"
        end if
        call nc_get(ncid, 'jlist', jlist)
        call nc_get(ncid, 'ixm_b', ixm)
        call nc_get(ncid, 'ixn_b', ixn)
        call nc_get(ncid, 'iota_b', iota_full)
        call nc_get(ncid, 'buco_b', buco_full)
        call nc_get(ncid, 'bvco_b', bvco_full)
        call nc_get(ncid, 'phi_b', phi_full)
        call nc_get(ncid, 'bmnc_b', bmnc_h)
        call nc_get(ncid, 'rmnc_b', rmnc_h)
        call nc_close(ncid)

        allocate (s_half(nsurf), rho_half(nsurf))
        do k = 1, nsurf
            s_half(k) = (real(jlist(k), dp) - 1.5_dp)/real(ns - 1, dp)
        end do
        rho_half = sqrt(s_half)

        torflux_si = -phi_full(ns)/TWOPI

        rho_min = 1.0e-3_dp
        rho_max = 1.0_dp
        s_min = rho_min**2
        s_max = rho_max**2

        d%n_rho = nrho
        d%n_s = nrho
        d%nfp = nfp_int
        d%torflux = torflux_si*GAUSS_CM2_PER_TM2
        d%rho_min = rho_min
        d%rho_max = rho_max
        d%h_rho = (rho_max - rho_min)/real(nrho - 1, dp)
        d%h_s = (s_max - s_min)/real(nrho - 1, dp)

        allocate (d%rho(nrho), d%s(nrho))
        allocate (d%A_phi(nrho), d%B_theta(nrho), d%B_phi(nrho))
        do ir = 1, nrho
            d%rho(ir) = rho_min + d%h_rho*real(ir - 1, dp)
            d%s(ir) = s_min + d%h_s*real(ir - 1, dp)
        end do

        call interp_1d(buco_full(jlist), rho_half, d%rho, nsurf, nrho, d%B_theta)
        call interp_1d(bvco_full(jlist), rho_half, d%rho, nsurf, nrho, d%B_phi)
        call iota_integral(iota_full(jlist), rho_half, d%s, nsurf, nrho, &
                           torflux_si, d%A_phi)
        d%B_theta = d%B_theta*GAUSS_CM_PER_TM
        d%B_phi = d%B_phi*GAUSS_CM_PER_TM
        d%A_phi = d%A_phi*GAUSS_CM2_PER_TM2

        allocate (theta(ntheta), zeta(nzeta), bmnc_out(nrho, nmn))
        do it = 1, ntheta
            theta(it) = TWOPI*real(it - 1, dp)/real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = TWOPI/real(nfp_int, dp)*real(iz - 1, dp)/real(nzeta, dp)
        end do

        call interp_modes(bmnc_h, rho_half, d%rho, ixm, nmn, nsurf, nrho, bmnc_out)

        allocate (bmod_geom(nrho, ntheta, nzeta))
        do ir = 1, nrho
            do it = 1, ntheta
                do iz = 1, nzeta
                    bmod_geom(ir, it, iz) = 0.0_dp
                    do mn = 1, nmn
                        angle = real(ixm(mn), dp)*theta(it) &
                                - real(ixn(mn), dp)*zeta(iz)
                        bmod_geom(ir, it, iz) = bmod_geom(ir, it, iz) &
                                                + bmnc_out(ir, mn)*cos(angle)
                    end do
                    bmod_geom(ir, it, iz) = bmod_geom(ir, it, iz)*GAUSS_PER_T
                end do
            end do
        end do

        d%n_theta = ntheta + 1
        d%n_phi = nzeta + 1
        d%h_theta = theta(2) - theta(1)
        d%h_phi = zeta(2) - zeta(1)
        allocate (d%Bmod(nrho, d%n_theta, d%n_phi))
        d%Bmod(:, 1:ntheta, 1:nzeta) = bmod_geom
        d%Bmod(:, d%n_theta, 1:nzeta) = bmod_geom(:, 1, :)
        d%Bmod(:, 1:ntheta, d%n_phi) = bmod_geom(:, :, 1)
        d%Bmod(:, d%n_theta, d%n_phi) = bmod_geom(:, 1, 1)

        mn00 = 0
        do mn = 1, nmn
            if (ixm(mn) == 0 .and. ixn(mn) == 0) then
                mn00 = mn
                exit
            end if
        end do
        if (mn00 > 0) then
            rmajor_m = rmnc_h(mn00, nsurf)
        else
            rmajor_m = 1.0_dp
        end if
        d%rmajor = rmajor_m

        call build_boozer_from_chartmap(d)
    end subroutine load_boozer_from_boozmn

    !> Interpolate boozmn Fourier coefficients from the half grid to rho_out.
    !> bmnc_h is shaped (nmn, nsurf): Fortran order from NetCDF (comput_surfs, mn_mode).
    subroutine interp_modes(bmnc_h, rho_half, rho_out, ixm, nmn, nsurf, nrho_out, &
                            bmnc_out)
        integer, intent(in) :: nmn, nsurf, nrho_out
        real(dp), intent(in) :: bmnc_h(nmn, nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: rho_out(nrho_out)
        integer, intent(in) :: ixm(nmn)
        real(dp), intent(out) :: bmnc_out(nrho_out, nmn)

        integer :: ir, mn, k
        real(dp) :: rho, frac, ratio

        do ir = 1, nrho_out
            rho = rho_out(ir)
            if (rho <= rho_half(1)) then
                do mn = 1, nmn
                    if (rho_half(1) > 0.0_dp) then
                        ratio = rho/rho_half(1)
                        bmnc_out(ir, mn) = bmnc_h(mn, 1)*ratio**min(ixm(mn), 50)
                    else
                        bmnc_out(ir, mn) = bmnc_h(mn, 1)
                    end if
                end do
            else if (rho >= rho_half(nsurf)) then
                do mn = 1, nmn
                    bmnc_out(ir, mn) = bmnc_h(mn, nsurf)
                end do
            else
                k = 1
                do while (k < nsurf)
                    if (rho_half(k + 1) >= rho) exit
                    k = k + 1
                end do
                if (rho_half(k + 1) > rho_half(k)) then
                    frac = (rho - rho_half(k))/(rho_half(k + 1) - rho_half(k))
                else
                    frac = 0.0_dp
                end if
                do mn = 1, nmn
                    bmnc_out(ir, mn) = (1.0_dp - frac)*bmnc_h(mn, k) &
                                       + frac*bmnc_h(mn, k + 1)
                end do
            end if
        end do
    end subroutine interp_modes

    !> Interpolate a 1D radial profile from the half grid to rho_out.
    subroutine interp_1d(vals_h, rho_half, rho_out, nsurf, nrho_out, vals_out)
        integer, intent(in) :: nsurf, nrho_out
        real(dp), intent(in) :: vals_h(nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: rho_out(nrho_out)
        real(dp), intent(out) :: vals_out(nrho_out)

        integer :: ir, k
        real(dp) :: rho, frac

        do ir = 1, nrho_out
            rho = rho_out(ir)
            if (rho <= rho_half(1)) then
                vals_out(ir) = vals_h(1)
            else if (rho >= rho_half(nsurf)) then
                vals_out(ir) = vals_h(nsurf)
            else
                k = 1
                do while (k < nsurf)
                    if (rho_half(k + 1) >= rho) exit
                    k = k + 1
                end do
                if (rho_half(k + 1) > rho_half(k)) then
                    frac = (rho - rho_half(k))/(rho_half(k + 1) - rho_half(k))
                else
                    frac = 0.0_dp
                end if
                vals_out(ir) = (1.0_dp - frac)*vals_h(k) + frac*vals_h(k + 1)
            end if
        end do
    end subroutine interp_1d

    !> A_phi(s) = -torflux_si * integral_0^s iota(s') ds', in T*m^2.
    subroutine iota_integral(iota_h, rho_half, s_out, nsurf, nrho_out, torflux_si, &
                             A_phi_out)
        integer, intent(in) :: nsurf, nrho_out
        real(dp), intent(in) :: iota_h(nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: s_out(nrho_out)
        real(dp), intent(in) :: torflux_si
        real(dp), intent(out) :: A_phi_out(nrho_out)

        integer :: ir, k
        real(dp) :: iota_at_s, s, s_half_sq(nsurf), cum(nsurf)

        s_half_sq = rho_half**2

        cum(1) = 0.0_dp
        do k = 2, nsurf
            cum(k) = cum(k - 1) + 0.5_dp*(iota_h(k - 1) + iota_h(k)) &
                     *(s_half_sq(k) - s_half_sq(k - 1))
        end do

        do ir = 1, nrho_out
            s = s_out(ir)
            if (s <= s_half_sq(1)) then
                iota_at_s = iota_h(1)
                A_phi_out(ir) = -torflux_si*iota_at_s*s
            else if (s >= s_half_sq(nsurf)) then
                A_phi_out(ir) = -torflux_si*(cum(nsurf) + &
                                             iota_h(nsurf)*(s - s_half_sq(nsurf)))
            else
                k = 1
                do while (k < nsurf)
                    if (s_half_sq(k + 1) >= s) exit
                    k = k + 1
                end do
                if (s_half_sq(k + 1) > s_half_sq(k)) then
                    iota_at_s = iota_h(k) + (iota_h(k + 1) - iota_h(k))* &
                                (s - s_half_sq(k))/(s_half_sq(k + 1) - s_half_sq(k))
                else
                    iota_at_s = iota_h(k)
                end if
                A_phi_out(ir) = -torflux_si*(cum(k) + 0.5_dp*(iota_h(k) + iota_at_s)* &
                                             (s - s_half_sq(k)))
            end if
        end do
    end subroutine iota_integral

end module boozmn_reader
