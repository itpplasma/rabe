program test_precession_analytic
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use mock_field_3d, only: mock_field_3d_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use precession, only: compute_precession_correction
    use utils, only: not_same, linspace

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.001_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field

    integer, parameter :: n_fieldlines = 6
    real(dp), parameter :: phi_tol = 5e-6

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: nfp = N_tor
    real(dp) :: iota
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp), parameter :: l_c = 1e-4, Omega_hat = 1e-2, s_tor = 0.25_dp
    real(dp), parameter :: expected_correction = 1.0_dp
    real(dp), parameter :: tol = 1e-10
    real(dp) :: correction
    logical :: test_failed

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

    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    call compute_precession_correction(field, fieldlines, &
                                       l_c, Omega_hat, s_tor, &
                                       correction)

    if (not_same(correction, expected_correction, abstol_in=tol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_precession_analytic failed: correction"
        print *, "found: ", correction
        print *, "expected: ", expected_correction
        test_failed = .true.
    end if

    if (test_failed) error stop

end program test_precession_analytic
