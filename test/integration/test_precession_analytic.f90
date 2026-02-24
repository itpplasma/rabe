program test_precession_analytic
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use precession, only: compute_precession_correction
    use utils, only: not_same, linspace

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field
    real(dp), parameter :: B_pert = 0.001_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field

    integer, parameter :: n_fieldlines = 50
    real(dp), parameter :: phi_tol = 1e-6

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.00_dp ! analytic formula for small iota
    real(dp), parameter :: nfp = 1.0_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp), parameter :: l_c = 1e-4, Omega_hat = 1e-2
    real(dp), parameter :: expected_correction = 1.0_dp
    real(dp), parameter :: tol = 1e-10
    real(dp) :: correction
    logical :: test_failed

    test_failed = .false.

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field%mock_perturbed_field_init(field, &
                                                   M_pol_pert, &
                                                   N_tor_pert, &
                                                   B_pert)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  iota, &
                                  perturbed_field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    call compute_precession_correction(field, l_c, Omega_hat, correction)

    if (not_same(correction, expected_correction, abstol_in=tol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_precession_analytic failed: correction"
        print *, "found: ", correction
        print *, "expected: ", expected_correction
        test_failed = .true.
    end if

    if (test_failed) error stop

end program test_precession_analytic
