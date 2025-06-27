module plot_quantities
    use constants, only: dp, pi
    use myplot_module, only: myplot
    use fieldline_mod, only: fieldline_t
    use utils, only: linspace

    implicit none

    type external_data_t
        real(dp), dimension(:), allocatable :: x
        real(dp), dimension(:), allocatable :: y
        character(len=1024) :: label
    end type external_data_t

contains

    subroutine plot_maxima_over_label(fieldlines)
        type(fieldline_t), dimension(:) :: fieldlines

        type(myplot) :: plt
        real(dp), dimension(size(fieldlines)) :: theta_0_shifted

        call plt%initialize(xlabel="$\vartheta_0 - \iota_p$ [$\pi$]", &
                            ylabel="$B$ [T]", &
                            figsize=(/7, 7/), &
                            legend=.true.)
        theta_0_shifted = modulo(fieldlines%theta_0 - fieldlines%iota_p, 2.0_dp*pi)
        call plt%add_plot(theta_0_shifted/pi, &
                          fieldlines%B_max(1), &
                          label="left $B_{max}$", &
                          linestyle="-o")
        call plt%add_plot(theta_0_shifted/pi, &
                          fieldlines%B_max(2), &
                          label="right $B_{max}$", &
                          linestyle="-o")
        call plt%show()
    end subroutine plot_maxima_over_label

    subroutine plot_field_along_chi_line(field, chi, M_pol, N_tor)
        use field_base, only: field_t
        use utils, only: linspace
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: chi, M_pol, N_tor

        type(myplot) :: plt
        integer, parameter :: n_points = 100
        real(dp), dimension(:), allocatable :: theta, phi, B

        character(len=1024) :: label
        integer :: current

        allocate (theta(n_points), phi(n_points), B(n_points))

        call linspace(0.0_dp, 2.0_dp*pi, n_points, theta)
        phi = (theta*M_pol - chi)/N_tor

        do current = 1, n_points
            call field%compute_B_mod(theta(current), phi(current), B(current))
        end do

        call plt%initialize(xlabel="$\vartheta$ [$\pi$]", &
                            ylabel="$B$ [T]", &
                            figsize=(/7, 7/), &
                            legend=.true.)
        write (label, "(A17,ES10.3E2)") "$B$ along $\chi=$", chi
        call plt%add_plot(theta/pi, B, label=label, linestyle="b-")
        call plt%show()

        deallocate (theta, phi, B)
    end subroutine plot_field_along_chi_line

    subroutine plot_fieldlines_over_field(fieldlines, field, N_tor)
        use field_base, only: field_t

        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: N_tor

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

    end subroutine plot_fieldlines_over_field

    subroutine plot_delta_eta(fieldlines, delta_eta_1, iota_p)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in), optional :: delta_eta_1, iota_p

        real(dp), dimension(size(fieldlines)) :: theta_0
        type(myplot) :: plt

        theta_0 = fieldlines%theta_0

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$\Delta \eta$", &
                            legend=.true.)

        call plt%add_plot(theta_0, &
                          fieldlines%delta_eta, &
                          label="$\Delta \eta$", &
                          linestyle="-")

        if (present(delta_eta_1)) then
            if (.not. present(iota_p)) error stop
            call plt%add_plot(theta_0, &
                              abs(delta_eta_1) - delta_eta_1*cos(theta_0 - iota_p), &
                              label="$\Delta \eta$ approx analytic", &
                              linestyle="-")
        end if

        call plt%show()
    end subroutine plot_delta_eta

    subroutine plot_delta_A(fieldlines, delta_A_1_analytic)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in), optional :: delta_A_1_analytic

        type(myplot) :: plt

        call plt%initialize(xlabel="$\vartheta_{mid} [\pi]$", &
                            ylabel="$\Delta A$", &
                            legend=.true.)

        call plt%add_plot(fieldlines%theta_0/pi, &
                          fieldlines%delta_aspect_ratio, &
                          label="$\Delta A$", &
                          linestyle="-")
        if (present(delta_A_1_analytic)) then
            call plt%add_plot(fieldlines%theta_0/pi, &
                              delta_A_1_analytic*cos(fieldlines%theta_0), &
                              label="$\Delta A$ approx analytic", &
                              linestyle="-")
        end if
        call plt%show()
    end subroutine plot_delta_A

    subroutine plot_I(fieldlines, I_0_analytic, I_1_analytic, eps_0, eps_1)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: I_0_analytic, I_1_analytic, eps_0, eps_1

        type(myplot) :: plt
        integer :: n_fieldlines
        real(dp) :: eps_ratio, I_factor, I_mean
        real(dp), dimension(size(fieldlines)) :: theta_0
        character(len=100) :: label

        n_fieldlines = size(fieldlines)
        theta_0 = fieldlines%theta_0
        eps_ratio = eps_1/abs(eps_0)
        I_mean = sum(fieldlines%integral_lambda_b_over_B_squared)/n_fieldlines
        I_factor = sum(fieldlines%integral_lambda_b_over_B_squared/ &
                       sqrt(1.0_dp - eps_ratio*cos(theta_0)))/n_fieldlines

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$I$ [T$^{-2}$]", &
                            legend=.true.)
        call plt%add_plot(theta_0, &
                          fieldlines%integral_lambda_b_over_B_squared, &
                          label="$I_{numeric}$", &
                          linestyle="b-")
        write (label, "(A23,F10.8)") "$I^{mean}_{numeric} = $", I_mean
        call plt%add_plot(theta_0, &
                          I_mean*cos(0.0_dp*theta_0), &
                          label=label, &
                          linestyle="b--")
        call plt%add_plot(theta_0, &
                          I_0_analytic + I_1_analytic*cos(theta_0), &
                          label="approx $I_{analytic}$", &
                          linestyle="r-")
        call plt%add_plot(theta_0, &
                          I_0_analytic*cos(0.0_dp*theta_0), &
                          label="approx $I^{mean}_{analytic}$", &
                          linestyle="r--")
        call plt%add_plot(theta_0, &
                          I_0_analytic*sqrt(1.0_dp - eps_ratio*cos(theta_0)), &
                          label="$I_{analytic}$", &
                          linestyle="g--")
        call plt%show()

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$I_{normalized}$ [1]", &
                            legend=.true.)
        call plt%add_plot(theta_0, &
                        fieldlines%integral_lambda_b_over_B_squared/I_factor - 1.0_dp, &
                          label="$I_{numeric}$", &
                          linestyle="b-")
        call plt%add_plot(theta_0, &
                          sqrt(1.0_dp - eps_ratio*cos(theta_0)) - 1.0_dp, &
                          label="$I_{analytic}$", &
                          linestyle="g--")
        call plt%show()
    end subroutine plot_I

    subroutine plot_I_factor(fieldlines, eps_0, eps_1)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: eps_0, eps_1

        type(myplot) :: plt
        real(dp), dimension(size(fieldlines)) :: theta_0, I_factor
        integer :: n_fieldlines
        real(dp) :: I_factor_mean
        character(len=1024) :: label

        theta_0 = fieldlines%theta_0
        I_factor = fieldlines%integral_lambda_b_over_B_squared &
                   /sqrt(1.0_dp - eps_1/abs(eps_0)*cos(theta_0))
        n_fieldlines = size(fieldlines)
        I_factor_mean = sum(I_factor)/n_fieldlines

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$I/\sqrt{1+\epsilon\cos{\vartheta_0}}$ [1/T$^2$]", &
                            legend=.true.)
        write (label, "(A38,F10.8)") "deviation from numeric $I_{factor} =$ ", &
            I_factor_mean
        call plt%add_plot(theta_0, &
                          I_factor/I_factor_mean, &
                          label=label, &
                          linestyle="k-")
        label = "analytic approximation $+ \mathcal{O}(\epsilon_0^2)$"
        call plt%add_plot(theta_0, &
                          -8.0_dp/(3.0_dp - 2.0_dp*eps_0)*eps_1*cos(theta_0) &
                          + 1.0_dp, &
                          label=label, &
                          linestyle="g-")
        call plt%show()
    end subroutine plot_I_factor

    subroutine plot_deviation(off_factor_A, off_factor_B, &
                              off_factor_A_analytic, off_factor_B_analytic, &
                              lambda_off_external)
        real(dp), intent(in) :: off_factor_A, off_factor_B

        real(dp), intent(in), optional :: off_factor_A_analytic
        real(dp), intent(in), optional :: off_factor_B_analytic

        type(external_data_t), intent(in), optional :: lambda_off_external

        type(myplot) :: plt
        integer, parameter :: n_points = 20
        real(dp), dimension(:), allocatable :: nu_star
        real(dp), dimension(2) :: nu_star_lim

        character(len=1024) :: label

        call plt%initialize(xlabel="$\nu_*$", &
                            ylabel="$|\lambda_{bB}|$", &
                            legend=.true., &
                            figsize=(/10, 10/))

        if (present(lambda_off_external)) then
            call plt%add_plot(lambda_off_external%x, &
                              lambda_off_external%y, &
                              label=lambda_off_external%label, &
                              linestyle="co", &
                              markersize=8, &
                              xscale="log", &
                              yscale="log")
            allocate (nu_star(size(lambda_off_external%x)), &
                      source=lambda_off_external%x)
        else
            allocate (nu_star(n_points))
            call linspace(0.0_dp, 8.0_dp, n_points, nu_star)
            nu_star = 0.1_dp**nu_star
        end if

        write (label, "(A35,ES10.3E2,A1)") "offset due to $\Delta A$ (factor$=$", &
            off_factor_A, ")"
        call plt%add_plot(nu_star, &
                          abs(off_factor_A)/sqrt(nu_star), &
                          label=label, &
                          linestyle="r-", &
                          xscale="log", &
                          yscale="log")

        if (present(off_factor_A_analytic)) then
            write (label, "(A42,ES10.3E2,A1)") &
                "analytic estimate (relative difference =", &
                abs(off_factor_A_analytic/off_factor_A - 1.0_dp), ")"
            call plt%add_plot(nu_star, &
                              abs(off_factor_A_analytic)/sqrt(nu_star), &
                              label=label, &
                              linestyle="ro", &
                              markersize=8, &
                              xscale="log", &
                              yscale="log")
        end if

        write (label, "(A38,ES10.3E2,A1)") "offset due to $\Delta \eta$ (factor$=$", &
            off_factor_B, ")"
        call plt%add_plot(nu_star, &
                          abs(off_factor_B)/nu_star, &
                          label=label, &
                          linestyle="b-", &
                          xscale="log", &
                          yscale="log")

        if (present(off_factor_B_analytic)) then
            write (label, "(A42,ES10.3E2,A1)") &
                "analytic estimate (relative difference =", &
                abs(off_factor_B_analytic/off_factor_B - 1.0_dp), ")"
            call plt%add_plot(nu_star, &
                              abs(off_factor_B_analytic)/nu_star, &
                              label=label, &
                              linestyle="bo", &
                              markersize=8, &
                              xscale="log", &
                              yscale="log")
        end if

        nu_star_lim(1) = minval(nu_star)
        nu_star_lim(2) = maxval(nu_star)
        call plt%add_plot(nu_star, &
                          abs(off_factor_A/sqrt(nu_star) + off_factor_B/nu_star), &
                          label="total", &
                          linestyle="c--", &
                          xscale="log", &
                          yscale="log", &
                          xlim=nu_star_lim)

        call plt%show()
    end subroutine plot_deviation

end module plot_quantities
