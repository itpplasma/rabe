program test_plot_anti_sigma
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_mod, only: make_flock_of_fieldlines

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = 0.0125_dp, eps_1 = 0.0005_dp
    type(anti_sigma_field_t) :: field
    real(dp), parameter :: B_pert = 0.01_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field

    real(dp), parameter :: phi_tol = 1e-4
    integer, parameter :: n_fieldlines = 20

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 1.0_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    call field%anti_sigma_field_init(N_tor, B_0, eps_0, eps_1)
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

    call plot_fieldlines_over_field(fieldlines, perturbed_field)

contains

    subroutine plot_fieldlines_over_field(fieldlines, field)
        use myplot_module, only: myplot
        use field_base, only: field_t

        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        class(field_t), intent(in) :: field

        type(myplot) :: plt
        integer :: current
        real(dp), dimension(2) :: phi_max, theta_max
        real(dp), dimension(2) :: phi_limits
        integer, parameter :: n_points = 100
        real(dp), dimension(n_points) :: phi, theta
        real(dp), dimension(n_points, n_points) :: phi_mesh, theta_mesh
        real(dp), dimension(n_points, n_points) :: B_mesh
        integer :: theta_idx, phi_idx
        character(len=100) :: label

        call plt%initialize(xlabel="$\varphi$", ylabel="$\vartheta$")

        phi_limits = fieldlines(1)%phi_0 + (/-1.5_dp*pi, +1.5_dp*pi/)
        do current = 1, size(fieldlines)
            phi_max = fieldlines(current)%phi_max
            call linspace(phi_max(1), &
                          phi_max(2), &
                          n_points, &
                          phi)
            theta = fieldlines(current)%get_theta(phi)
            write (label, '(F0.1)') fieldlines(current)%theta_0
            call plt%add_plot(phi, theta, &
                              label=label, &
                              linestyle="k.", &
                              linewidth=1)
            theta_max = fieldlines(current)%get_theta(phi_max)
            call plt%add_plot(phi_max, theta_max, &
                              label=label, &
                              linestyle="ro", &
                              markersize=5)
        end do

        call linspace(phi_limits(1), phi_limits(2), n_points, phi)
        call linspace(0.0_dp, 2.0_dp*pi, n_points, theta)
        do theta_idx = 1, n_points
            do phi_idx = 1, n_points
                call field%compute_B_mod(theta(theta_idx), &
                                         phi(phi_idx), &
                                         B_mesh(theta_idx, phi_idx))
            end do
        end do
        call plt%add_contour(phi, theta, B_mesh, &
                             levels=20, &
                             colorbar=.true., &
                             filled=.true.)
        call plt%show()

    end subroutine plot_fieldlines_over_field

end program test_plot_anti_sigma
