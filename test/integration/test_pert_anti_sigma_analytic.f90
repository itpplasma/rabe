program test_pert_anti_sigma_analytic
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use fieldline_integrals, only: fourier_transform_over_label
    use fieldline_integrals, only: fieldline_modes_t
    use utils, only: linspace
    use utils, only: not_same

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = 0.0125_dp, eps_1 = 0.0005_dp
    type(anti_sigma_field_t) :: field
    real(dp), parameter :: B_pert = 0.001_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field
    real(dp), parameter :: B_max = B_0*(1.0_dp + eps_0)
    real(dp), parameter :: delta_eta_1 = -B_pert/B_max**2.0_dp
    real(dp), parameter :: delta_eta_0 = B_pert/B_max**2.0_dp

    real(dp), parameter :: reltol_delta_eta = (B_pert/B_max)
    real(dp), parameter :: abstol = 1e-12
    real(dp), parameter :: phi_tol = 2e-5

    integer, parameter :: n_fieldlines = 50
    integer, parameter :: n_modes = n_fieldlines/2 + 1

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 1.0_dp ! analytic formula for small iota
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines
    type(fieldline_modes_t) :: fieldline_modes

    integer :: current
    real(dp) :: B_mod
    real(dp), dimension(:), allocatable :: zeros
    real(dp) :: expected_deviation

    logical :: test_failed

    test_failed = .false.

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field%mock_perturbed_field_init(field, &
                                                   M_pol_pert, &
                                                   N_tor_pert, &
                                                   B_pert)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, &
                                  theta_0, &
                                  iota, &
                                  perturbed_field, &
                                  M_pol, &
                                  N_tor, &
                                  phi_tol)
    call fourier_transform_over_label(fieldlines, &
                                      fieldline_modes)

    if (not_same(delta_eta_1, fieldline_modes%delta_eta%cos_coeffs(2), &
                 reltol_in=reltol_delta_eta, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_pert_anti_sigma_analytic failed: 1st delta_eta cos mode"
        print *, "found: ", fieldline_modes%delta_eta%cos_coeffs(2)
        print *, "expected: ", delta_eta_1
        print *, "ratio: ", fieldline_modes%delta_eta%cos_coeffs(2)/delta_eta_1
        test_failed = .true.
    end if

    if (not_same(delta_eta_0, fieldline_modes%delta_eta%cos_coeffs(1), &
                 reltol_in=reltol_delta_eta)) then
        print *, "-------------------------------------------------------------"
        print *, "test_pert_anti_sigma_analytic failed: 0th delta_eta cos mode"
        print *, "found: ", fieldline_modes%delta_eta%cos_coeffs(1)
        print *, "expected: ", delta_eta_0
        print *, "ratio: ", fieldline_modes%delta_eta%cos_coeffs(1)/delta_eta_0
        test_failed = .true.
    end if

    allocate (zeros(size(fieldline_modes%radial_drift%cos_coeffs)))
    zeros = 0.0_dp

    if (not_same(zeros, fieldline_modes%delta_eta%sin_coeffs, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_pert_anti_sigma_analytic failed: delta_eta sin modes"
        print *, "found: ", fieldline_modes%delta_eta%sin_coeffs
        print *, "expected: all ", zeros(1)
        test_failed = .true.
    end if

    deallocate (zeros)

    do current = 3, n_modes
        expected_deviation = B_pert/B_max*(0.5_dp*B_pert/B_max)**(current - 2)
        if (not_same(0.0_dp, &
                     fieldline_modes%delta_eta%cos_coeffs(current), &
                     abstol_in=max(abstol, expected_deviation) &
                     )) then
            print *, "-------------------------------------------------------------"
            print *, "test_pert_anti_sigma_analytic failed: ", &
                "delta_eta cos mode number", current - 1
            print *, "found: ", fieldline_modes%delta_eta%cos_coeffs(current)
            print *, "expected: ", expected_deviation, "or lower"
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_pert_anti_sigma_analytic
