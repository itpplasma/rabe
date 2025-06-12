program test_deviation_helical_anti_sigma
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation

    implicit none

    real(dp), parameter :: M_pol = 2.0_dp, N_tor = 10.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.005, eps_1 = -0.000375_dp
    real(dp), parameter :: factor = sqrt(abs(eps_0)/(1.0_dp + abs(eps_0))) &
                           /B_0**2.0_dp/N_tor*sqrt(32.0_dp)
    real(dp), parameter :: I_0_analytic = factor*(1.0_dp + abs(eps_0)*2.0_dp/3.0_dp)
    real(dp), parameter :: eps_ratio = eps_1/abs(eps_0)
    real(dp), parameter :: I_1_analytic = factor*(-0.5_dp*eps_ratio)* &
                           (1.0_dp + 6.0_dp*abs(eps_0))
    real(dp), parameter :: B_max = B_0*(1.0_dp + abs(eps_0))
    real(dp), parameter :: delta_A_1 = -0.25_dp*abs(eps_1/eps_0)
    real(dp), parameter :: psi_edge = abs(-0.28274_dp)/(2.0_dp*pi) !Tm^2
    real(dp), parameter :: dr_dpsi = 1.0_dp/0.0661777_dp/psi_edge/100.0_dp
    real(dp), parameter :: J_pol_over_N_tor = -4.0_dp*1e6
    real(dp), parameter :: I_tor = 0.0_dp
    real(dp), parameter :: R = 8.0_dp
    type(anti_sigma_field_t) :: field
    real(dp), dimension(2), parameter :: chi_max = (/+pi, -pi/) + 2.0_dp*pi

    real(dp), parameter :: phi_tol = 1e-6
    real(dp), parameter :: B_max_reltol = phi_tol**2.0_dp*N_tor**2.0
    real(dp), parameter :: I_mean_reltol = eps_0*eps_0*3.0_dp
    real(dp), parameter :: I_amplitude_reltol = max(eps_ratio*eps_ratio*5.0_dp, &
                                                    eps_0*eps_0*10.0_dp)
    integer, parameter :: n_fieldlines = 20
    real(dp), dimension(2) :: phi_max

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.0_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    integer :: current
    real(dp), dimension(n_fieldlines) :: I_shifted
    real(dp) :: I_mean, I_amplitude
    logical :: test_failed

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)
    call make_flock_of_fieldlines(fieldlines, &
                                  theta_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  phi_tol)

    test_failed = .false.
    do current = 1, n_fieldlines
        phi_max = (M_pol*fieldlines(current)%theta_0 - chi_max)/N_tor
        if (not_same(fieldlines(current)%phi_max(1), &
                     phi_max(1), &
                     abstol_in=2.0_dp*phi_tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_deviation_helical_anti_sigma failed: left phi_max"
            print *, "found: ", fieldlines(current)%phi_max(1)
            print *, "analytic: ", phi_max(1)
            test_failed = .true.
        end if
        if (not_same(fieldlines(current)%phi_max(2), &
                     phi_max(2), &
                     abstol_in=2.0_dp*phi_tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_deviation_helical_anti_sigma failed: right phi_max"
            print *, "found: ", fieldlines(current)%phi_max(2)
            print *, "analytic: ", phi_max(2)
            test_failed = .true.
        end if
     if (not_same(1.0_dp/fieldlines(current)%eta_b, B_max, reltol_in=B_max_reltol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_deviation_helical_anti_sigma failed: B_max"
            print *, "found: ", 1.0_dp/fieldlines(current)%eta_b
            print *, "analytic: ", B_max
            print *, "relative error: ", fieldlines(current)%eta_b*B_max - 1.0_dp
            test_failed = .true.
        end if
    end do

    I_mean = sum(fieldlines%integral_lambda_b_over_B_squared)/n_fieldlines
    if (not_same(I_mean, I_0_analytic, reltol_in=I_mean_reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_deviation_helical_anti_sigma failed: I_mean"
        print *, "found: ", I_mean
        print *, "analytic: ", I_0_analytic
        print *, "relative error: ", I_mean/I_0_analytic - 1.0_dp
        test_failed = .true.
    end if

    I_shifted = fieldlines%integral_lambda_b_over_B_squared - I_mean
    I_amplitude = 0.5_dp*(maxval(I_shifted) - minval(I_shifted))

    if (not_same(I_amplitude, I_1_analytic, &
                 reltol_in=I_amplitude_reltol, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_deviation_helical_anti_sigma failed: I_amplitude"
        print *, "found: ", I_amplitude
        print *, "analytic: ", I_1_analytic
        print *, "relative error: ", I_amplitude/I_1_analytic - 1.0_dp
        test_failed = .true.
    end if

    if (test_failed) error stop

    !call plot_fieldlines_over_field(fieldlines, field)
    !call plot_delta_eta(fieldlines)
    call plot_I_factor(fieldlines, eps_0, eps_1)
    !call plot_I(fieldlines, I_0_analytic, I_1_analytic, eps_0, eps_1)
    !call plot_delta_A(fieldlines, delta_A_1)

    !call calc_deviation(fieldlines, field, deviation_A, deviation_B)

    !covariant_factor = -2.0_dp*1e-7*(J_pol_over_N_tor*abs(N_tor) + I_tor*iota)
    !off_factor_A = deviation_A*dr_dpsi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    !off_factor_B = deviation_B*0.5*R*pi*dr_dpsi

    !call plot_deviation(off_factor_A, off_factor_B)

contains

    subroutine plot_delta_eta(fieldlines)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines

        type(myplot) :: plt

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$\Delta \eta$")

        call plt%add_plot(fieldlines%theta_0, &
                          fieldlines%delta_eta, &
                          label="\Delta \eta", &
                          linestyle="-")
        call plt%show()
    end subroutine plot_delta_eta

    subroutine plot_delta_A(fieldlines, delta_A_1_analytic)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: delta_A_1_analytic

        type(myplot) :: plt

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$\Delta A$", &
                            legend=.true.)

        call plt%add_plot(fieldlines%theta_0, &
                          fieldlines%delta_aspect_ratio, &
                          label="$\Delta A$", &
                          linestyle="-")
        call plt%add_plot(fieldlines%theta_0, &
                          delta_A_1_analytic*cos(fieldlines%theta_0), &
                          label="$\Delta A$ analytic", &
                          linestyle="-")
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
                            ylabel="$I/\sqrt{1 + \epsilon\cos{\vartheta_0}}$ [T$^{-2}$]", &
                            legend=.true.)
        write (label, "(A38,F10.8)") "deviation from numeric $I_{factor} =$ ",// &
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

    subroutine plot_deviation(off_factor_A, off_factor_B)
        real(dp), intent(in) :: off_factor_A, off_factor_B

        type(myplot) :: plt
        integer, parameter :: n_points = 100
        real(dp), dimension(n_points) :: nu_star

        call plt%initialize(xlabel="$\nu_*$", &
                            ylabel="$\lambda_{bB}$", &
                            legend=.true.)
        call linspace(0.0_dp, 8.0_dp, n_points, nu_star)
        nu_star = 0.1_dp**nu_star

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
    end subroutine plot_deviation

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

    end subroutine plot_fieldlines_over_field

end program test_deviation_helical_anti_sigma
