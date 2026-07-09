!> \file ad_driver.f90
!! \brief Flat, NetCDF-free `bind(C)` entry point into the FourierField bootstrap
!! pipeline, for the flang/Enzyme AD backend (`rabe.ad`).
!!
!! \details This mirrors the per-surface workflow of `app/main.f90` /
!! `python/example_fourier.py`, but exposes it as a single C-callable procedure over
!! plain arrays and scalars — no Fortran derived types cross the boundary. That keeps
!! the flang extension free of f90wrap/opaque handles (so it cannot be confused with
!! the gfortran extension objects) and matches the shape Enzyme differentiates: the
!! eventual adjoint is just another `bind(C)` symbol over the same buffers.
module ad_driver
    use iso_c_binding, only: c_int, c_double
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
    !! \details Builds a `fourier_field_t` from the modes, constructs the flock of
    !! field lines, and evaluates the offset coefficients and critical collisionality
    !! — the same numbers `python/example_fourier.py` obtains via
    !! `FlockOfFieldlines.calc_offset_coefficients` / `calc_nu_star_crit`.
    !!
    !! \param[in]  mn_max            number of Fourier modes
    !! \param[in]  m                 poloidal mode numbers (length mn_max)
    !! \param[in]  n                 toroidal mode numbers, normalised to nfp (length mn_max)
    !! \param[in]  B_mn              Fourier coefficients of B in Tesla (length mn_max)
    !! \param[in]  B_theta_cov       covariant poloidal component of B in T*m
    !! \param[in]  B_phi_cov         covariant toroidal component of B in T*m
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
        mn_max, m, n, B_mn, B_theta_cov, B_phi_cov, nfp, n_grid, &
        iota, M_pol, N_tor, max_n_fieldlines, R, dr_dAtheta, &
        lambda_a, lambda_b, nu_star_crit) &
        bind(C, name="rabe_fourier_offset_coefficients")
        integer(c_int), value :: mn_max
        integer(c_int), intent(in) :: m(mn_max), n(mn_max)
        real(c_double), intent(in) :: B_mn(mn_max)
        real(c_double), value :: B_theta_cov, B_phi_cov
        integer(c_int), value :: nfp, n_grid, max_n_fieldlines
        real(c_double), value :: iota, M_pol, N_tor, R, dr_dAtheta
        real(c_double), intent(out) :: lambda_a, lambda_b, nu_star_crit

        type(fourier_field_t) :: field
        type(flock_of_fieldlines_t) :: flock

        call fourier_field_init(field, m, n, B_mn, B_theta_cov, B_phi_cov, &
                                nfp, n_grid)
        call make_flock_of_fieldlines(flock, max_n_fieldlines, iota, field, &
                                      M_pol, N_tor, real(nfp, dp))
        call calc_offset_coefficients(flock, R, dr_dAtheta, lambda_a, lambda_b)
        nu_star_crit = calc_nu_star_crit(flock, R)

    end subroutine rabe_fourier_offset_coefficients

end module ad_driver
