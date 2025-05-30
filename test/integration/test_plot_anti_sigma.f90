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
    real(dp), parameter :: phi_0 = pi
    type(anti_sigma_field_t) :: field
    real(dp), parameter :: B_pert = 0.1_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field

    real(dp), parameter :: phi_tol = 1e-4
    integer, parameter :: n_fieldlines = 50

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.0_dp
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
        real(dp), dimension(2) :: phi, theta
        character(len=100) :: label

        call plt%initialize()
        do current = 1, size(fieldlines)
            theta = fieldlines(current)%theta_0
            phi = fieldlines(current)%phi_max
            write (label, '(F0.1)') fieldlines(current)%theta_0
            call plt%add_plot(phi, theta, label=label, linestyle="-")
        end do
        call plt%show()

    end subroutine plot_fieldlines_over_field

end program test_plot_anti_sigma
