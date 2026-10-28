program test_zero_precession
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use mock_field_3d, only: mock_field_3d_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use precession, only: compute_precession_correction
    use precession, only: get_smallest_maximum
    use precession, only: get_biggest_minimum
    use precession, only: set_fieldline_minima
    use precession_reference, only: integrate_radial_drift
    use grid_mod, only: integration_grid_t
    use grid_mod, only: fieldline_for_precession_t
    use grid_mod, only: set_integration_grids
    use grid_mod, only: compute_bounce_integrals
    use grid_mod, only: set_splines
    use fieldline_labels, only: fourier_transform_over_label
    use fieldline_labels, only: fieldline_modes_t
    use fourier, only: real_ft
    use utils, only: not_same, linspace

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, &
                           eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.001_dp, M_pol_pert = 1.0_dp, &
                           N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field

    integer, parameter :: n_fieldlines = 6

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: nfp = N_tor
    real(dp) :: iota
    type(flock_of_fieldlines_t) :: flock

    real(dp), parameter :: l_c = 1e-4, s_tor = 0.25_dp
    real(dp), parameter :: Omega_hat_zero = 0.0_dp
    real(dp), parameter :: tol = 1e-6
    real(dp) :: correction
    real(dp), dimension(:), allocatable :: flux_modes
    type(fieldline_modes_t) :: fieldline_modes
    integer :: n_modes, idx
    logical :: test_failed

    real(dp) :: B_covariant_phi, B_covariant_theta
    real(dp) :: conversion_factor
    real(dp) :: expected, found
    real(dp) :: dummy(11), hcovar(3), hctrvr(3)

    ! numerical check variables
    type(fieldline_for_precession_t), dimension(:), allocatable :: &
        prec_fieldlines
    type(integration_grid_t) :: grid
    real(dp) :: eta_t, eta_c
    real(dp), dimension(n_fieldlines) :: radial_drift_integrals
    integer, parameter :: n_modes_ft = n_fieldlines/2 + 1
    real(dp), dimension(n_modes_ft) :: numerical_cos, numerical_sin

    test_failed = .false.

    call field_2D%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field_2D%mock_perturbed_field_init(field_2D, &
                                                      M_pol_pert, &
                                                      N_tor_pert, &
                                                      B_pert)
    call field%mock_field_3d_init(perturbed_field_2D)
    call field%get_iota(s_tor, iota)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call make_flock_from_labels(flock, &
                                xi_0, &
                                iota, &
                                field, &
                                M_pol, &
                                N_tor, &
                                nfp)

    call fourier_transform_over_label(flock, fieldline_modes)

    call field%evaluate([s_tor, 0.0_dp, 0.0_dp], dummy(1), dummy(2), &
                        dummy(3:5), hcovar, hctrvr, dummy(9:11))
    B_covariant_phi = hcovar(3)
    B_covariant_theta = hcovar(2)
    conversion_factor = (B_covariant_phi + iota*B_covariant_theta) &
                        /(0.75_dp*field%psi_tor_edge)

    call compute_precession_correction(field, flock, &
                                       l_c, Omega_hat_zero, s_tor, &
                                       correction, flux_modes)

    n_modes = n_fieldlines/2 + 1

    ! --- Test 1: compare flux_modes with fieldline integral modes ---
    do idx = 2, n_modes
        expected = fieldline_modes%radial_drift%sin_coeffs(idx) &
                   *conversion_factor
        found = flux_modes(idx)
        if (not_same(found, expected, reltol_in=tol, &
                     abstol_in=0.0_dp)) then
            print *, "---------------------------------------------------"
            print *, "test_precession_analytic failed: flux_mode", idx
            print *, "found: ", found
            print *, "expected: ", expected
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_zero_precession
