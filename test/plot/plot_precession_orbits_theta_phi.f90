program plot_precession_orbits_theta_phi
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
    integer, parameter :: n_orbits_plot = 6
    integer, dimension(n_orbits_plot), parameter :: selected = [1, 5, 10, 20, 30, 40]
    real(dp), parameter :: gauss2tesla = 1e-4_dp

    type(myplot) :: plt
    type(myplot) :: plt_orbit
    integer :: i, idx, n_steps, step, unit_id
    real(dp), dimension(:), allocatable :: phi, theta
    real(dp), dimension(:), allocatable :: lambda_orbit, lambda_sq_orbit
    real(dp), dimension(:), allocatable :: B_orbit, B_orbit_tesla
    real(dp), dimension(:), allocatable :: inv_eta_line, invariant_check
    real(dp), dimension(:), allocatable :: step_axis
    real(dp) :: bmod, sqrtg
    real(dp), dimension(3) :: x, bder, hcovar, hctrvr, hcurl
    character(len=128) :: label
    character(len=128) :: orbit_file
    character(len=128) :: orbit_title
    logical :: no_bounce

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
                        title="Selected trace_orbit_till_bounce paths")

    do i = 1, n_orbits_plot
        idx = selected(i)
        if (idx < 1 .or. idx > size(orbit_history, 3)) cycle
        n_steps = 0
        do step = 1, size(orbit_history, 1)
            if (ieee_is_nan(orbit_history(step, 1, idx))) exit
            n_steps = step
        end do
        if (n_steps < 2) cycle

        allocate (phi(n_steps), theta(n_steps), lambda_orbit(n_steps), lambda_sq_orbit(n_steps), &
                  B_orbit(n_steps), B_orbit_tesla(n_steps), inv_eta_line(n_steps), &
                  invariant_check(n_steps), step_axis(n_steps))
        phi = orbit_history(1:n_steps, 3, idx)
        theta = orbit_history(1:n_steps, 2, idx)
        lambda_orbit = orbit_history(1:n_steps, 5, idx)
        lambda_sq_orbit = lambda_orbit**2.0_dp

        do step = 1, n_steps
            x = orbit_history(step, 1:3, idx)
            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            B_orbit(step) = bmod
            B_orbit_tesla(step) = bmod*gauss2tesla
            inv_eta_line(step) = 1.0_dp/upper_grid%eta(idx)
            invariant_check(step) = B_orbit_tesla(step)*upper_grid%eta(idx) + &
            lambda_sq_orbit(step)
            step_axis(step) = real(step, dp)
        end do

        no_bounce = ieee_is_nan(upper_grid%I_j(idx))
        if (no_bounce) then
      write (label, '(A,I0,A,F8.5)') 'idx=', idx, ' no-bounce eta=', upper_grid%eta(idx)
        else
        write (label, '(A,I0,A,F8.5)') 'idx=', idx, ' bounced eta=', upper_grid%eta(idx)
        end if

        write (orbit_file, '(A,I0,A)') 'orbit_theta_phi_idx_', idx, '.dat'
        open (newunit=unit_id, file=trim(orbit_file), status='replace', action='write')
        write (unit_id, '(A)') '# step s theta phi p lambda I_j_integral tau_integral'
        do step = 1, n_steps
         write (unit_id, '(I8,1X,7(ES23.15E3,1X))') step, orbit_history(step, 1, idx), &
orbit_history(step, 2, idx), orbit_history(step, 3, idx), orbit_history(step, 4, idx), &
   orbit_history(step, 5, idx), orbit_history(step, 6, idx), orbit_history(step, 7, idx)
        end do
        close (unit_id)

        write (orbit_file, '(A,I0,A)') 'orbit_pitch_bmod_idx_', idx, '.dat'
        open (newunit=unit_id, file=trim(orbit_file), status='replace', action='write')
        write (unit_id, '(A)') '# step lambda lambda_sq bmod_tesla inv_eta bmod_tesla*eta+lambda_sq'
        do step = 1, n_steps
            write (unit_id, '(I8,1X,5(ES23.15E3,1X))') step, lambda_orbit(step), lambda_sq_orbit(step), &
                B_orbit_tesla(step), inv_eta_line(step), invariant_check(step)
        end do
        close (unit_id)

    write (orbit_title, '(A,I0,A,F8.5)') 'Orbit idx=', idx, ' eta=', upper_grid%eta(idx)
        call plt_orbit%initialize(xlabel='step', &
                                  ylabel='value', &
                                  legend=.true., &
                                  title=trim(orbit_title))
        call plt_orbit%add_plot(step_axis, (1.0_dp - lambda_orbit**2)/inv_eta_line, label='(1-lambda^2)/eta', linestyle='b.')
     call plt_orbit%add_plot(step_axis, B_orbit_tesla, label='Bmod [T]', linestyle='r.')
     call plt_orbit%add_plot(step_axis, inv_eta_line, label='1/eta [T]', linestyle='k.')
        call plt_orbit%add_plot(step_axis, invariant_check, label='B*eta + lambda^2', linestyle='g.')
        call plt_orbit%show()

        call plt%add_plot(phi, theta, label=trim(label), linestyle='.')
        deallocate (phi, theta, lambda_orbit, lambda_sq_orbit, B_orbit, B_orbit_tesla, &
                    inv_eta_line, invariant_check, step_axis)
    end do

    call plt%show()

end program plot_precession_orbits_theta_phi
