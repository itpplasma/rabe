program test_lambda_helical_anti_sigma
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use fieldline_integrands, only: lambda_over_B_squared

    implicit none

    real(dp), parameter :: M_pol = 2.0_dp, N_tor = 10.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.05, eps_1 = -0.00375_dp
    type(anti_sigma_field_t) :: field

    real(dp), parameter :: phi_tol = 1e-6
    integer, parameter :: n_fieldlines = 20

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.0_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp), dimension(n_fieldlines) :: lambda_integral
    integer :: current
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

    lambda_integral = sqrt(abs(eps_0)/(1.0_dp + abs(eps_0))) &
                      /N_tor*4.0_dp*sqrt(2.0_dp) &
                      *sqrt(1.0_dp - eps_1/abs(eps_0)*cos(theta_0))

    test_failed = .false.
    do current = 1, n_fieldlines
        if (.true.) then
            print *, "-------------------------------------------------------------"
            print *, "test_lambda_helical_anti_sigma failed: lambda"
            print *, "found: "
            print *, "analytic: "
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

contains

    subroutine plot_I(fieldlines, I_0_analytic, I_1_analytic, eps_0, eps_1)
        use myplot_module, only: myplot
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: I_0_analytic, I_1_analytic, eps_0, eps_1

        integer :: n_fieldlines
        type(myplot) :: plt
        real(dp) :: I_mean
        real(dp), dimension(size(fieldlines)) :: theta_0

        n_fieldlines = size(fieldlines)
        theta_0 = fieldlines%theta_0
        I_mean = sum(fieldlines%integral_lambda_b_over_B_squared)/n_fieldlines

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$I$ [T$^{-2}$]", &
                            legend=.true.)
        call plt%add_plot(theta_0, &
                          fieldlines%integral_lambda_b_over_B_squared, &
                          label="$I_{numeric}$", &
                          linestyle="b-")
        call plt%add_plot(theta_0, &
                          I_mean*cos(0.0_dp*theta_0), &
                          label="$I^{mean}_{numeric}$", &
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
                          I_0_analytic*sqrt(1.0_dp - eps_1/abs(eps_0)*cos(theta_0)), &
                          label="$I_{analytic}$", &
                          linestyle="g--")
        call plt%show()

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$I_{normalized}$ [1]", &
                            legend=.true.)
        call plt%add_plot(theta_0, &
                          fieldlines%integral_lambda_b_over_B_squared/I_mean - 1.0_dp, &
                          label="$I_{numeric}$", &
                          linestyle="b-")
        call plt%add_plot(theta_0, &
                          sqrt(1.0_dp - eps_1/abs(eps_0)*cos(theta_0)) - 1.0_dp, &
                          label="$I_{analytic}$", &
                          linestyle="g--")
        call plt%show()
    end subroutine plot_I

end program test_lambda_helical_anti_sigma
