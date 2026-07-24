!> \file ad_driver.f90
!! \brief Single entry point wrapping the FourierField bootstrap pipeline, from a
!! flat list of Boozer Fourier modes to the asymptotic offset coefficients.
!!
!! \details Mirrors the per-surface workflow of `app/main.f90` /
!! `python/example_fourier.py` (`fourier_field_init` -> `make_flock_of_fieldlines` ->
!! `calc_offset_coefficients`/`calc_nu_star_crit`), collapsed into one procedure so
!! there is a single, well-defined target for the eventual Enzyme adjoint. Wrapped
!! for Python via f90wrap alongside the rest of `rabe_lib` (see
!! `python/CMakeLists.txt`) - no special calling convention needed; Enzyme
!! differentiates the compiled procedure directly regardless of how it is exposed
!! to Python.
module ad_driver

    use constants, only: dp
    use fourier_field, only: fourier_field_t, fourier_field_init
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines
    use coefficients, only: calc_offset_coefficients, calc_nu_star_crit

    implicit none

contains

    !>
    !! \brief Compute the asymptotic bootstrap offset coefficients for a magnetic
    !! field given as a flat list of Boozer Fourier modes.
    !!
    !! \param[in]  m                 poloidal mode numbers (flat array, length mn_max)
    !! \param[in]  n                 toroidal mode numbers, normalised to nfp (length mn_max)
    !! \param[in]  B_mn              Fourier coefficients of B in Tesla (length mn_max)
    !! \param[in]  B_theta_covariant covariant poloidal component of B in T*m
    !! \param[in]  B_phi_covariant   covariant toroidal component of B in T*m
    !! \param[in]  nfp               number of field periods
    !! \param[in]  n_grid            spline grid points per angle direction
    !! \param[in]  iota              rotational transform of the surface
    !! \param[in]  M_pol             dominant poloidal helicity mode number
    !! \param[in]  N_tor             dominant toroidal helicity mode number (== nfp or 0)
    !! \param[in]  max_n_fieldlines  upper bound on number of field lines
    !! \param[in]  R                 major radius in metres
    !! \param[in]  dr_dAtheta        dr/dA_theta conversion factor [rad/(T m)]
    !! \param[out] lambda_a          1/sqrt(nu_star) offset coefficient
    !! \param[out] lambda_b          1/nu_star offset coefficient
    !! \param[out] nu_star_crit      lower collisionality validity limit
    !<
    subroutine rabe_fourier_offset_coefficients( &
        m, n, B_mn, B_theta_covariant, B_phi_covariant, nfp, n_grid, &
        iota, M_pol, N_tor, max_n_fieldlines, R, dr_dAtheta, &
        lambda_a, lambda_b, nu_star_crit)
        integer, intent(in) :: m(:), n(:)
        real(dp), intent(in) :: B_mn(:)
        real(dp), intent(in) :: B_theta_covariant, B_phi_covariant
        integer, intent(in) :: nfp
        integer, intent(in) :: n_grid
        real(dp), intent(in) :: iota
        real(dp), intent(in) :: M_pol, N_tor
        integer, intent(in) :: max_n_fieldlines
        real(dp), intent(in) :: R, dr_dAtheta
        real(dp), intent(out) :: lambda_a, lambda_b, nu_star_crit

        type(fourier_field_t) :: field
        type(flock_of_fieldlines_t) :: flock

        call fourier_field_init(field, m, n, B_mn, B_theta_covariant, B_phi_covariant, &
                                nfp, n_grid)
        call make_flock_of_fieldlines(flock, max_n_fieldlines, iota, field, &
                                      M_pol, N_tor, real(nfp, dp))
        call calc_offset_coefficients(flock, R, dr_dAtheta, lambda_a, lambda_b)
        nu_star_crit = calc_nu_star_crit(flock, R)

    end subroutine rabe_fourier_offset_coefficients

end module ad_driver
