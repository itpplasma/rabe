module synthetic_spectra
    !! Small synthetic Boozer spectra for the derivative-mode tests: a
    !! quasi-helical (QH), a quasi-axisymmetric (QA) and a quasi-isodynamic
    !! (QI) surface with a few symmetry-breaking modes, so that all rabe
    !! outputs are nonzero.
    use constants, only: dp
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t

    implicit none
    private

    public :: synthetic_qh, synthetic_qa, synthetic_qi, shift_half_period

contains

    subroutine synthetic_qh(spec, settings)
        type(surface_spectrum_t), intent(out) :: spec
        type(derivative_settings_t), intent(out) :: settings

        ! columns: m, n/nfp, bmnc, nabla_s_mnc
        real(dp), parameter :: table(4, 9) = reshape([ &
                                                     0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, &
                                                     1.0_dp, 1.0_dp, -0.1_dp, 0.3_dp, &
                                                     0.0_dp, 1.0_dp, 4.0e-3_dp, 0.05_dp, &
                                                     1.0_dp, 0.0_dp, 3.0e-3_dp, 0.1_dp, &
                                                     2.0_dp, 2.0_dp, -8.0e-3_dp, 0.02_dp, &
                                                     2.0_dp, 1.0_dp, 2.0e-3_dp, 0.0_dp, &
                                                     1.0_dp, 2.0_dp, 1.5e-3_dp, 0.0_dp, &
                                                     1.0_dp, -1.0_dp, 1.0e-3_dp, 0.01_dp, &
                                                     3.0_dp, 2.0_dp, 5.0e-4_dp, 0.0_dp], [4, 9])

        ! Helicity M_pol = -1, N_tor = nfp as in the golden QH case, so the
        ! dominant mode is cos(theta + nfp*phi): flip the sign of n.
        call fill(spec, table, nfp=4.0_dp)
        spec%xn = -spec%xn
        spec%iota = 1.137_dp
        spec%G = 10.0_dp
        spec%I = 0.12_dp
        spec%psi_edge = 1.5_dp
        spec%R = 10.0_dp
        settings%M_pol = -1.0_dp
        settings%N_tor = 4.0_dp
        settings%max_n_fieldlines = 60
        settings%n_gauss = 64
        settings%sign_sqrtg = 1
    end subroutine synthetic_qh

    subroutine synthetic_qa(spec, settings)
        type(surface_spectrum_t), intent(out) :: spec
        type(derivative_settings_t), intent(out) :: settings

        real(dp), parameter :: table(4, 8) = reshape([ &
                                                     0.0_dp, 0.0_dp, 1.0_dp, 1.5_dp, &
                                                     1.0_dp, 0.0_dp, -0.12_dp, 0.25_dp, &
                                                     0.0_dp, 1.0_dp, 5.0e-3_dp, 0.04_dp, &
                                                     1.0_dp, 1.0_dp, 4.0e-3_dp, 0.02_dp, &
                                                     1.0_dp, -1.0_dp, -2.0e-3_dp, 0.0_dp, &
                                                     2.0_dp, 0.0_dp, 6.0e-3_dp, 0.01_dp, &
                                                     2.0_dp, 1.0_dp, 1.0e-3_dp, 0.0_dp, &
                                                     0.0_dp, 2.0_dp, 1.0e-3_dp, 0.0_dp], [4, 8])

        call fill(spec, table, nfp=2.0_dp)
        spec%iota = 0.4217_dp
        spec%G = 8.0_dp
        spec%I = -0.05_dp
        spec%psi_edge = -0.9_dp
        spec%R = 6.0_dp
        settings%M_pol = 1.0_dp
        settings%N_tor = 0.0_dp
        settings%max_n_fieldlines = 60
        settings%n_gauss = 64
        settings%sign_sqrtg = -1
    end subroutine synthetic_qa

    !> QI-like surface (M_pol = 0, N_tor = nfp) in booz_xform convention: the
    !! dominant mirror term cos(nfp*phi) puts the maximum of B at the origin
    !! and the well centre at phi = pi/nfp, as in SQuID and Nautilus wouts.
    subroutine synthetic_qi(spec, settings)
        type(surface_spectrum_t), intent(out) :: spec
        type(derivative_settings_t), intent(out) :: settings

        real(dp), parameter :: table(4, 9) = reshape([ &
                                                     0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, &
                                                     0.0_dp, 1.0_dp, 0.15_dp, 0.05_dp, &
                                                     1.0_dp, 0.0_dp, 2.0e-2_dp, 0.3_dp, &
                                                     1.0_dp, 1.0_dp, 3.0e-2_dp, 0.02_dp, &
                                                     1.0_dp, -1.0_dp, -1.0e-2_dp, 0.01_dp, &
                                                     0.0_dp, 2.0_dp, 2.0e-2_dp, 0.0_dp, &
                                                     2.0_dp, 1.0_dp, 3.0e-3_dp, 0.0_dp, &
                                                     1.0_dp, 2.0_dp, 2.0e-3_dp, 0.0_dp, &
                                                     2.0_dp, 0.0_dp, -2.0e-3_dp, 0.01_dp], [4, 9])

        call fill(spec, table, nfp=4.0_dp)
        spec%iota = 0.8731_dp
        spec%G = 17.0_dp
        spec%I = 0.05_dp
        spec%psi_edge = -1.2_dp
        spec%R = 17.0_dp
        settings%M_pol = 0.0_dp
        settings%N_tor = 4.0_dp
        settings%max_n_fieldlines = 60
        settings%n_gauss = 64
        settings%sign_sqrtg = -1
    end subroutine synthetic_qi

    !> Move the toroidal origin by half a field period, phi -> phi + pi/nfp:
    !! the coefficients of modes with odd xn/nfp change sign.
    subroutine shift_half_period(spec, sign_k)
        type(surface_spectrum_t), intent(inout) :: spec
        real(dp), allocatable, intent(out) :: sign_k(:)

        sign_k = 1.0_dp - 2.0_dp*modulo(nint(spec%xn/spec%nfp), 2)
        spec%bmnc = sign_k*spec%bmnc
        spec%nabla_s_mnc = sign_k*spec%nabla_s_mnc
    end subroutine shift_half_period

    subroutine fill(spec, table, nfp)
        type(surface_spectrum_t), intent(inout) :: spec
        real(dp), intent(in) :: table(:, :), nfp

        spec%nfp = nfp
        spec%xm = table(1, :)
        spec%xn = table(2, :)*nfp
        spec%bmnc = table(3, :)
        spec%nabla_s_mnc = table(4, :)
    end subroutine fill

end module synthetic_spectra
