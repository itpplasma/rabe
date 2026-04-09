program plot_bounce_orbits
    use constants, only: dp, pi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
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
    use field_instance, only: initialize_field_instance, magfie
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
    real(dp), parameter :: phi_tol = 8e-6_dp
    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota
    real(dp), parameter :: nfp = 1.0_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines
    type(fieldline_with_minimum_t) :: precession_fieldline
    type(integration_grid_t) :: lower_grid, upper_grid

    real(dp), parameter :: s_tor = 0.25_dp
    real(dp) :: B_min, phi_bottom
    real(dp) :: lowest_B_max, eta_t, eta_c

    real(dp), dimension(:, :, :), allocatable :: orbit_history
    integer, parameter :: n_orbits_plot = 3
    integer, dimension(n_orbits_plot), parameter :: selected = [1, 5, 10]
    real(dp), parameter :: gauss2tesla = 1e-4_dp

    type(myplot) :: plt
    type(myplot) :: plt_orbit
    integer :: i, idx, n_steps, step
    real(dp), dimension(:), allocatable :: phi, theta
    real(dp), dimension(:), allocatable :: lambda
    real(dp), dimension(:), allocatable :: B_orbit
    real(dp), dimension(:), allocatable :: check
    real(dp), dimension(:), allocatable :: step_axis
    real(dp) :: eta_inv
    real(dp) :: bmod, sqrtg
    real(dp), dimension(3) :: x, bder, hcovar, hctrvr, hcurl
    character(len=128) :: label
    character(len=128) :: orbit_title
    integer, dimension(2), parameter :: figsize = [14, 10]

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

    call set_integration_grids(eta_t, eta_c, lower_grid, upper_grid)
    call initialize_field_instance(field)

    ! Capture orbit points for each eta level of the upper grid.
    call compute_bounce_integrals(field, precession_fieldline, s_tor, upper_grid, &
                                  orbit_history=orbit_history)

    call plt%initialize(xlabel="$\varphi$", &
                        ylabel="$\vartheta$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$\vartheta(\varphi)$ bounce-orbit traces")

    do i = 1, n_orbits_plot
        idx = selected(i)
        n_steps = 0
        do step = 1, size(orbit_history, 1)
            if (ieee_is_nan(orbit_history(step, 1, idx))) exit
            n_steps = step
        end do
        if (n_steps < 2) cycle

        allocate (phi(n_steps), theta(n_steps), lambda(n_steps), &
                  B_orbit(n_steps), check(n_steps), step_axis(n_steps))
        phi = orbit_history(1:n_steps, 3, idx)
        theta = orbit_history(1:n_steps, 2, idx)
        lambda = orbit_history(1:n_steps, 5, idx)
        eta_inv = 1.0_dp/upper_grid%eta(idx)

        do step = 1, n_steps
            x = orbit_history(step, 1:3, idx)
            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            B_orbit(step) = bmod*gauss2tesla
            check(step) = B_orbit(step)*upper_grid%eta(idx) + lambda(step)**2.0_dp
            step_axis(step) = real(step, dp)
        end do

        write (label, '(A,I0,A,F8.5,A)') '$i=', idx, &
            '$, $\eta=', upper_grid%eta(idx), '$'

        write (orbit_title, '(A,I0,A,F8.5,A)') '$i=', idx, &
            '$, $\eta=', upper_grid%eta(idx), '$'
        call plt_orbit%initialize(xlabel='$n$', &
                                  ylabel='$\mathrm{diagnostic}$', &
                                  legend=.true., &
                                  figsize=figsize, &
                                  title=trim(orbit_title))
        call plt_orbit%add_plot(step_axis, (1.0_dp - lambda**2)*eta_inv, &
                                label='$(1-\lambda^2)/\eta$', linestyle='bo')
        call plt_orbit%add_plot(step_axis, B_orbit, &
                                label='$B\,[\mathrm{T}]$', linestyle='r.')
        call plt_orbit%add_plot(step_axis, 0.0_dp*step_axis + eta_inv, &
                                label='$1/\eta\,[\mathrm{T}]$', linestyle='k.')
        call plt_orbit%add_plot(step_axis, check, &
                                label='$B\eta+\lambda^2$', linestyle='g.')
        call plt_orbit%show()

        call plt%add_plot(phi, theta, label=trim(label), linestyle='.')
        deallocate (phi, theta, lambda, B_orbit, check, step_axis)
    end do

    call plt%show()

end program plot_bounce_orbits
