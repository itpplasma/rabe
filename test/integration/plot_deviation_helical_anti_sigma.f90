program plot_deviation_helical_anti_sigma
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation

    implicit none

    real(dp), parameter :: M_pol = 2.0_dp, N_tor = 10.0_dp
    character(len=*), parameter :: bc_filename = "input/helical_anti.bc"
    real(dp), parameter :: psi_edge = abs(-0.28274_dp)/(2.0_dp*pi) !Tm^2
    real(dp), parameter :: dr_dpsi = 1.0_dp/0.0661777_dp/psi_edge/100.0_dp
    real(dp), parameter :: J_pol_over_N_tor = -4.0_dp*1e6
    real(dp), parameter :: I_tor = 0.0_dp
    real(dp), parameter :: stor = 0.9999_dp, R = 8.0_dp
    real(dp), parameter :: eps_0 = -0.05, eps_1 = -0.0375
    real(dp), parameter :: delta_A_1 = -0.25_dp*abs(eps_1/eps_0)
    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 7e-7
    integer, parameter :: n_fieldlines = 50

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.47_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp) :: deviation_A, deviation_B
    integer, parameter :: n_points = 100
    real(dp), dimension(n_points) :: nu_star
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B
    type(myplot) :: plt

    call field%neo_field_init(bc_filename, stor)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, &
                                  theta_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  phi_tol)

    call plot_fieldlines_over_field(fieldlines, field)

    call calc_deviation(fieldlines, field, deviation_A, deviation_B)

    call plt%initialize(xlabel="$\vartheta_{mid}$", &
                        ylabel="$\Delta \eta$")

    call plt%add_plot(fieldlines%theta_0, &
                      fieldlines%delta_eta, &
                      label="\Delta \eta", &
                      linestyle="-")
    call plt%show()

    call plt%initialize(xlabel="$\vartheta_{mid}$", &
                        ylabel="$\Delta A$", &
                        legend=.true.)

    call plt%add_plot(fieldlines%theta_0, &
                      fieldlines%delta_aspect_ratio, &
                      label="$\Delta A$", &
                      linestyle="-")
    call plt%add_plot(fieldlines%theta_0, &
                      delta_A_1*cos(fieldlines%theta_0), &
                      label="$\Delta A$ analytic", &
                      linestyle="-")
    call plt%show()

    call plt%initialize(xlabel="$\nu_*$", &
                        ylabel="$\lambda_{bB}$", &
                        legend=.true.)
    call linspace(0.0_dp, 8.0_dp, n_points, nu_star)
    nu_star = 0.1_dp**nu_star

    covariant_factor = -2.0_dp*1e-7*(J_pol_over_N_tor*abs(N_tor) + I_tor*iota)
    off_factor_A = deviation_A*dr_dpsi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dpsi

    call plt%add_plot(nu_star, &
                      off_factor_A/sqrt(nu_star), &
                      label="offset factor due to aspect ratio =", &
                      linestyle="r-", &
                      xscale="log", &
                      yscale="log")

    call plt%add_plot(nu_star, &
                      off_factor_B/nu_star, &
                      label="due to misaligment", &
                      linestyle="b-", &
                      xscale="log", &
                      yscale="log")

    call plt%show()

contains

    subroutine plot_fieldlines_over_field(fieldlines, field)
        use myplot_module, only: myplot
        use field_base, only: field_t

        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        class(field_t), intent(in) :: field

        type(myplot) :: plt
        integer :: current
        integer :: n_fieldlines
        real(dp), dimension(2) :: phi_max, theta_max
        real(dp), dimension(2) :: phi_limits
        integer, parameter :: n_points = 300
        real(dp), dimension(n_points) :: phi, theta
        real(dp), dimension(n_points, n_points) :: phi_mesh, theta_mesh
        real(dp), dimension(n_points, n_points) :: B_mesh
        integer :: theta_idx, phi_idx
        character(len=100) :: label

        n_fieldlines = size(fieldlines)

        call plt%initialize(xlabel="$\varphi$", ylabel="$\vartheta$")

        phi_limits = +(/-1.5_dp*pi/N_tor + fieldlines(1)%phi_0, &
                        +1.5_dp*pi/N_tor + fieldlines(n_fieldlines)%phi_0/)
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
        call plt%add_contour(phi, theta, transpose(B_mesh), &
                             levels=20, &
                             colorbar=.true., &
                             filled=.true.)
        call plt%show()

!    open (unit=10, file="B_mesh.dat", status='replace', action='write', form='formatted')
!         do theta_idx = 1, n_points
!          write (10, '( *(F11.9,1X) )') (B_mesh(theta_idx, phi_idx), phi_idx=1, n_points)
!         end do

    end subroutine plot_fieldlines_over_field

end program plot_deviation_helical_anti_sigma
