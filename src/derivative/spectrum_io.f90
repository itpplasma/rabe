module spectrum_io
    !! Spectrum input for the derivative mode: projection of any `field_t` on a
    !! flux surface onto a truncated Boozer cosine spectrum, and reading a
    !! spectrum from NetCDF (layout documented in src/derivative/README.md).
    use constants, only: dp, pi
    use field_base, only: field_t
    use spectrum_derivative, only: surface_spectrum_t

    implicit none
    private

    public :: boozer_mode_set, spectralize_field, read_spectrum_file

contains

    !> booz_xform mode ordering: m = 0 with n = 0..nboz, then m = 1..mboz with
    !! n = -nboz..nboz; `xn = n*nfp`.
    subroutine boozer_mode_set(mboz, nboz, nfp, xm, xn)
        integer, intent(in) :: mboz, nboz
        real(dp), intent(in) :: nfp
        real(dp), allocatable, intent(out) :: xm(:), xn(:)

        integer :: m, n, k

        allocate (xm((nboz + 1) + mboz*(2*nboz + 1)))
        allocate (xn(size(xm)))
        k = 0
        do m = 0, mboz
            do n = merge(0, -nboz, m == 0), nboz
                k = k + 1
                xm(k) = real(m, dp)
                xn(k) = real(n, dp)*nfp
            end do
        end do
    end subroutine boozer_mode_set

    !> Project B and |grad s| of a surface-fixed field onto cos(m*theta - n*phi)
    !! by an exact discrete Fourier transform on a 4(mboz+1) x 4(nboz+1) grid
    !! per field period. `residual` is the maximum of |B - B_spectrum|/B_00 on
    !! the grid shifted by half a cell, and `asymmetry` is max|bmns|/B_00.
    subroutine spectralize_field(field, nfp, mboz, nboz, spec, residual, asymmetry)
        use spectral_surface, only: spectral_surface_t, spectral_surface_init
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: nfp
        integer, intent(in) :: mboz, nboz
        type(surface_spectrum_t), intent(inout) :: spec
            !! iota, G, I, psi_edge, R are set by the caller
        real(dp), intent(out) :: residual, asymmetry

        integer :: n_theta, n_phi, i, j, k
        real(dp), allocatable :: theta(:), phi(:), B(:, :), S(:, :), bmns(:)
        real(dp) :: arg, weight, B_I, G_cov, B_spec
        type(spectral_surface_t) :: check

        n_theta = 4*(mboz + 1)
        n_phi = 4*(nboz + 1)
        allocate (theta(n_theta), phi(n_phi), B(n_theta, n_phi), S(n_theta, n_phi))
        theta = [(2.0_dp*pi*real(i - 1, dp)/real(n_theta, dp), i=1, n_theta)]
        phi = [(2.0_dp*pi*real(j - 1, dp)/(real(n_phi, dp)*nfp), j=1, n_phi)]
        do j = 1, n_phi
            do i = 1, n_theta
                call field%compute_B_mod(theta(i), phi(j), B(i, j))
                call field%compute_nabla_s(theta(i), phi(j), S(i, j))
            end do
        end do

        spec%nfp = nfp
        call boozer_mode_set(mboz, nboz, nfp, spec%xm, spec%xn)
        allocate (spec%bmnc(size(spec%xm)), spec%nabla_s_mnc(size(spec%xm)))
        allocate (bmns(size(spec%xm)))
        spec%bmnc = 0.0_dp
        spec%nabla_s_mnc = 0.0_dp
        bmns = 0.0_dp
        do k = 1, size(spec%xm)
            weight = 2.0_dp/real(n_theta*n_phi, dp)
            if (spec%xm(k) == 0.0_dp .and. spec%xn(k) == 0.0_dp) &
                weight = 0.5_dp*weight
            do j = 1, n_phi
                do i = 1, n_theta
                    arg = spec%xm(k)*theta(i) - spec%xn(k)*phi(j)
                    spec%bmnc(k) = spec%bmnc(k) + weight*B(i, j)*cos(arg)
                    spec%nabla_s_mnc(k) = spec%nabla_s_mnc(k) + weight*S(i, j)*cos(arg)
                    bmns(k) = bmns(k) + weight*B(i, j)*sin(arg)
                end do
            end do
        end do
        asymmetry = maxval(abs(bmns))/abs(spec%bmnc(1))

        call field%get_covariant_components(B_I, G_cov)
        call spectral_surface_init(check, spec%xm, spec%xn, spec%bmnc, &
                                   spec%nabla_s_mnc, B_I, G_cov)
        residual = 0.0_dp
        do j = 1, n_phi
            do i = 1, n_theta
                call field%compute_B_mod(theta(i) + pi/real(n_theta, dp), &
                                         phi(j) + pi/(real(n_phi, dp)*nfp), &
                                         B(i, j))
                call check%compute_B_mod(theta(i) + pi/real(n_theta, dp), &
                                         phi(j) + pi/(real(n_phi, dp)*nfp), &
                                         B_spec)
                residual = max(residual, abs(B(i, j) - B_spec))
            end do
        end do
        residual = residual/abs(spec%bmnc(1))
    end subroutine spectralize_field

    !> Read all surfaces of a spectrum file. Variables (Fortran order):
    !! xm(mode), xn(mode), nfp, s_tor(surface), bmnc(mode, surface),
    !! nabla_s_mnc(mode, surface), iota, G, I, psi_edge, R (surface).
    subroutine read_spectrum_file(filename, s_tor, spectra)
        use netcdf_mod, only: netcdf_t
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: s_tor(:)
        type(surface_spectrum_t), allocatable, intent(out) :: spectra(:)

        type(netcdf_t) :: nc
        integer :: n_surf, n_mode, i
        real(dp) :: nfp
        real(dp), allocatable :: xm(:), xn(:), bmnc(:, :), smnc(:, :)
        real(dp), allocatable :: iota(:), G(:), I_cov(:), psi(:), R(:)

        call nc%open(filename)
        call nc%read_dim("surface", n_surf)
        call nc%read_dim("mode", n_mode)
        allocate (xm(n_mode), xn(n_mode), bmnc(n_mode, n_surf), smnc(n_mode, n_surf))
        allocate (s_tor(n_surf), iota(n_surf), G(n_surf), I_cov(n_surf), &
                  psi(n_surf), R(n_surf))
        call nc%read_real("nfp", nfp)
        call nc%read_real_1d("xm", xm)
        call nc%read_real_1d("xn", xn)
        call nc%read_real_1d("s_tor", s_tor)
        call nc%read_real_2d("bmnc", bmnc)
        call nc%read_real_2d("nabla_s_mnc", smnc)
        call nc%read_real_1d("iota", iota)
        call nc%read_real_1d("G", G)
        call nc%read_real_1d("I", I_cov)
        call nc%read_real_1d("psi_edge", psi)
        call nc%read_real_1d("R", R)
        call nc%close()

        allocate (spectra(n_surf))
        do i = 1, n_surf
            spectra(i)%nfp = nfp
            spectra(i)%xm = xm
            spectra(i)%xn = xn
            spectra(i)%bmnc = bmnc(:, i)
            spectra(i)%nabla_s_mnc = smnc(:, i)
            spectra(i)%iota = iota(i)
            spectra(i)%G = G(i)
            spectra(i)%I = I_cov(i)
            spectra(i)%psi_edge = psi(i)
            spectra(i)%R = R(i)
        end do
    end subroutine read_spectrum_file

end module spectrum_io
