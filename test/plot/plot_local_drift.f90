program plot_anti_sigma
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: fieldline_t

    use plot_quantities, only: plot_fieldline_over_local_drift

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp, nfp = N_tor
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field
    real(dp), parameter :: B_pert = -0.01_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field

    real(dp), parameter :: phi_tol = 8e-5

    real(dp), parameter :: iota = 1.0_dp
    real(dp), parameter :: eta = 0.5_dp
    type(fieldline_t) :: fieldline

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field%mock_perturbed_field_init(field, &
                                                   M_pol_pert, &
                                                   N_tor_pert, &
                                                   B_pert)
    fieldline%iota = iota
    fieldline%nfp = nfp
    fieldline%N_tor = N_tor
    fieldline%M_pol = M_pol
    fieldline%phi_0 = pi/nfp
    fieldline%theta_0 = 0.0_dp

    call plot_fieldline_over_local_drift(fieldline, perturbed_field, eta)

end program plot_anti_sigma
