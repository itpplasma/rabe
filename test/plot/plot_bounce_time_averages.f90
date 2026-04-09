program plot_bounce_time_averages
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
    use precession, only: set_splines
    use precession, only: evaluate_grid_splines
    use field_instance, only: initialize_field_instance
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
    type(integration_grid_t) :: grid

    real(dp), parameter :: s_tor = 0.25_dp
    integer, parameter :: n_spline_plot = 300
    real(dp), dimension(:), allocatable :: eta
    real(dp), dimension(:), allocatable :: bounce_time
    real(dp), dimension(:), allocatable :: bounce_time_deep
    real(dp), dimension(:), allocatable :: I_j, I_j_boundary
    real(dp), dimension(:), allocatable :: t_spline
    real(dp), dimension(:), allocatable :: eta_spline
    real(dp), dimension(:), allocatable :: I_j_spline_vals
    real(dp), dimension(:), allocatable :: bounce_coef_spline_vals
    real(dp) :: B_min, eta_b, phi_bottom
    real(dp) :: lowest_B_max, eta_t, eta_c, theta_min
    real(dp) :: well_depth
    real(dp) :: dummy, sqrtg
    real(dp), dimension(3) :: x, b_der, h_covar, h_ctrvr, h_curl
    integer, dimension(2), parameter :: figsize = [14, 10]

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

    call set_integration_grids(eta_t, eta_c, grid)
    call initialize_field_instance(field)
    call compute_bounce_integrals(field, precession_fieldline, s_tor, grid)
    call set_splines(grid)

    allocate (eta(grid%n_grid - 1))
    allocate (bounce_time(size(eta)))
    allocate (I_j(size(eta)))

    eta = grid%eta(2:grid%n_grid)

    bounce_time = grid%bounce_time(2:grid%n_grid)

    I_j = grid%I_j(2:grid%n_grid)

    theta_min = precession_fieldline%get_theta(phi_bottom)
    x = [s_tor, theta_min, phi_bottom]
    call field%evaluate(x, dummy, sqrtg, b_der, h_covar, h_ctrvr, h_curl)

    allocate (bounce_time_deep(size(eta)))
    allocate (I_j_boundary(size(eta)))
    allocate (t_spline(n_spline_plot))
    allocate (eta_spline(n_spline_plot))
    allocate (I_j_spline_vals(n_spline_plot))
    allocate (bounce_coef_spline_vals(n_spline_plot))

    well_depth = max(machine_eps, 1.0_dp - B_min*eta_b)
    bounce_time_deep = 4.0_dp*pi/abs(h_ctrvr(2))/sqrt(well_depth)
    I_j_boundary = 4.0_dp*sqrt(well_depth)/h_ctrvr(3)
    call linspace(grid%t(1), grid%t(grid%n_grid), n_spline_plot, t_spline)
    eta_spline = eta_c + t_spline**2.0_dp
    call evaluate_grid_splines(grid, t_spline, I_j_spline_vals, bounce_coef_spline_vals)

    call plt%initialize(xlabel="$1 - \eta/\eta_t$", &
                        ylabel="$\tau_b v_{\mathrm{th}}$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$\tau_b$ vs $1 - \eta/\eta_t$")

    call plt%add_plot(1.0_dp - eta/eta_t, bounce_time, &
                      label="$\tau_b$ (numerical)", &
                      linestyle="r-")
    call plt%add_plot(1.0_dp - eta/eta_t, bounce_time_deep, &
                      label="$\tau_b$ (deep-trapped estimate)", &
                      linestyle="k--")
    call plt%show()

    call plt%initialize(xlabel="$\eta/\eta_c - 1$", &
                        ylabel="$I_j$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$I_j$ vs $\eta/\eta_c - 1$")

    call plt%add_plot(eta/eta_c - 1.0_dp, I_j, &
                      label="$I_j$ (numerical)", &
                      linestyle="r-")
    call plt%add_plot(eta/eta_c - 1.0_dp, I_j_boundary, &
                      label="$I_j$ (trapped-passing boundary estimate)", &
                      linestyle="k--")
    call plt%add_plot(eta_spline/eta_c - 1.0_dp, I_j_spline_vals, &
                      label="$I_j$ (spline)", &
                      linestyle="b--")

    call plt%show()

    call plt%initialize(xlabel="$t$", &
                        ylabel="$C_{\mathrm{bounce}}$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$C_{\mathrm{bounce}}$ vs $t$")

    call plt%add_plot(grid%t, grid%bounce_coef, &
                      label="$C_{\mathrm{bounce}}$ (numerical)", &
                      linestyle="r-")
    call plt%add_plot(t_spline, bounce_coef_spline_vals, &
                      label="$C_{\mathrm{bounce}}$ (spline)", &
                      linestyle="b--")

    call plt%show()

end program plot_bounce_time_averages
