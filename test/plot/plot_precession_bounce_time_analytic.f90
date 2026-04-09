program plot_precession_bounce_time_analytic
    use constants, only: dp, pi, machine_eps
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use mock_field_3d, only: mock_field_3d_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use precession, only: integration_grid_t
    use precession, only: fieldline_with_minimum_t
    use precession, only: get_fieldline_at_global_maximum
    use precession, only: find_magnetic_well_bottom
    use precession, only: set_integration_grids
    use precession, only: compute_bounce_integrals
    use field_instance, only: initialize_field_instance
    use plot_quantities, only: plot_local_drift_over_fieldline
    use utils, only: linspace
    use myplot_module, only: myplot

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.001_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field

    integer, parameter :: n_fieldlines = 50
    real(dp), parameter :: phi_tol = 5e-6_dp
    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota
    real(dp), parameter :: nfp = N_tor
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines
    type(fieldline_with_minimum_t) :: precession_fieldline
    type(integration_grid_t) :: lower_grid, upper_grid

    real(dp), parameter :: s_tor = 0.25_dp
    real(dp), dimension(:), allocatable :: eta
    real(dp), dimension(:), allocatable :: eta_level_plot
    real(dp), dimension(:), allocatable :: bounce_time
    real(dp), dimension(:), allocatable :: bounce_time_deep
    real(dp), dimension(:), allocatable :: I_j, I_j_boundary
    real(dp) :: B_min, eta_b, phi_bottom
    real(dp) :: lowest_B_max, eta_t, eta_c, theta_min
    real(dp) :: well_depth
    real(dp) :: dummy, sqrtg
    real(dp), dimension(3) :: x, b_der, h_covar, h_ctrvr, h_curl

    type(myplot) :: plt

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

    precession_fieldline%fieldline_t = get_fieldline_at_global_maximum(fieldlines)
    call find_magnetic_well_bottom(field, precession_fieldline, phi_bottom, B_min)
    lowest_B_max = minval(precession_fieldline%B_max)
    eta_t = 1.0_dp/B_min
    eta_c = 1.0_dp/lowest_B_max
    precession_fieldline%phi_min = phi_bottom
    precession_fieldline%B_min = B_min
    eta_b = precession_fieldline%eta_b

    call set_integration_grids(eta_t, eta_c, lower_grid, upper_grid)
    call initialize_field_instance(field)
    call compute_bounce_integrals(field, precession_fieldline, s_tor, lower_grid)
    call compute_bounce_integrals(field, precession_fieldline, s_tor, upper_grid)

    eta = lower_grid%eta(2:lower_grid%n_grid)
    bounce_time = lower_grid%bounce_time(2:lower_grid%n_grid)
    I_j = lower_grid%I_j(2:lower_grid%n_grid)

    theta_min = precession_fieldline%get_theta(phi_bottom)
    x = [s_tor, theta_min, phi_bottom]
    call field%evaluate(x, dummy, sqrtg, b_der, h_covar, h_ctrvr, h_curl)

    allocate (bounce_time_deep(size(eta)))
    allocate (I_j_boundary(size(eta)))

    well_depth = max(machine_eps, 1.0_dp - B_min*eta_b)
    bounce_time_deep = 4.0_dp*pi/abs(h_ctrvr(2))/sqrt(well_depth)
    I_j_boundary = 4.0_dp*sqrt(well_depth)/h_ctrvr(3)

    call plt%initialize(xlabel="$1 - \eta/\eta_t$", &
                        ylabel="$\tau_b v_{th}$", &
                        legend=.true., &
                        title="Bounce time")

    call plt%add_plot(1.0_dp - eta/eta_t, bounce_time, &
                      label="numerical", &
                      linestyle="k-")
    call plt%add_plot(1.0_dp - eta/eta_t, bounce_time_deep, &
                      label="deep trapped estimate", &
                      linestyle="r--")
    call plt%show()

    call plt%initialize(xlabel="$\eta/\eta_c - 1$", &
                        ylabel="$I_j$", &
                        legend=.true., &
                        title="adiabatic invariant")

    call plt%add_plot(eta/eta_c - 1.0_dp, I_j, &
                      label="numerical", &
                      linestyle="k-")
    call plt%add_plot(eta/eta_c - 1.0_dp, I_j_boundary, &
                      label="trapped-passing boundary estimate", &
                      linestyle="r--")

    call plt%show()

    allocate (eta_level_plot((lower_grid%n_grid - 1) + (upper_grid%n_grid - 1)))
    eta_level_plot(1:lower_grid%n_grid - 1) = lower_grid%eta(2:lower_grid%n_grid)
    eta_level_plot(lower_grid%n_grid:) = upper_grid%eta(2:upper_grid%n_grid)

    call plot_local_drift_over_fieldline(field, precession_fieldline%fieldline_t, eta_level_plot, &
                                         interval=precession_fieldline%phi_max)

end program plot_precession_bounce_time_analytic
