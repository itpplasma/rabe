module synthetic_spectra
    !! Small synthetic Boozer spectra for the derivative-mode tests: a
    !! quasi-helical (QH) and a quasi-axisymmetric (QA) surface with a few
    !! symmetry-breaking modes, so that all rabe outputs are nonzero.
    use constants, only: dp
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t

    implicit none
    private

    public :: synthetic_qh, synthetic_qa

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
