program test_lambda_helical_anti_sigma
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use fieldline_labels, only: get_theta_0
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use integrate, only: integrate_1d_substituted

    implicit none

    real(dp), parameter :: M_pol = 2.0_dp, N_tor = 10.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.05, eps_1 = -0.00375_dp
    type(anti_sigma_field_t) :: field

    real(dp), parameter :: phi_tol = 1e-6
    integer, parameter :: max_n_fieldlines = 20

    real(dp), dimension(:), allocatable :: xi_0
    real(dp), parameter :: iota = 0.0_dp, nfp = max(1.0_dp, abs(N_tor))
    real(dp) :: approx_iota
    integer :: n_fieldlines
    type(fieldline_t), dimension(:), allocatable :: fieldlines
    type(fieldline_t) :: current_fieldline

    real(dp), dimension(:), allocatable :: lambda_integral_analytic
    real(dp), dimension(:), allocatable :: lambda_integral
    real(dp), parameter :: lambda_reltol = 1e-8
    integer :: current
    logical :: test_failed

    logical :: should_plot = .false.

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call get_theta_0(max_n_fieldlines, iota, M_pol, N_tor, nfp, xi_0, approx_iota)
    n_fieldlines = size(xi_0)
    allocate (fieldlines(n_fieldlines))
    allocate (lambda_integral_analytic(n_fieldlines))
    allocate (lambda_integral(n_fieldlines))
    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  approx_iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    lambda_integral_analytic = sqrt(abs(eps_0)/(1.0_dp + abs(eps_0))) &
                               /N_tor*4.0_dp*sqrt(2.0_dp) &
                               *sqrt(1.0_dp - eps_1/abs(eps_0)*cos(fieldlines%theta_0))

    lambda_integral = 0.0_dp
    do current = 1, n_fieldlines
        current_fieldline = fieldlines(current)
        call integrate_1d_substituted(wrapper_lambda, &
                                      current_fieldline%phi_max(1), &
                                      current_fieldline%phi_max(2), &
                                      lambda_integral(current))
    end do

    test_failed = .false.
    if (not_same(lambda_integral_analytic, &
                 lambda_integral, &
                 reltol_in=lambda_reltol, &
                 abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_lambda_helical_anti_sigma failed: lambda"
        print *, "found: ", lambda_integral
        print *, "analytic: ", lambda_integral_analytic
        test_failed = .true.
    end if

    if (should_plot) then
        call plot_lambda(xi_0, lambda_integral_analytic, lambda_integral)
    end if

    if (test_failed) error stop

contains

    function wrapper_lambda(phi)
        use fieldline_integrands, only: lambda_over_B_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_lambda

        real(dp) :: theta, B

        theta = current_fieldline%get_theta(phi)
        wrapper_lambda = lambda_over_B_squared(field, &
                                               theta, &
                                               phi, &
                                               current_fieldline%eta_b)
        call field%compute_B_mod(theta, phi, B)
        wrapper_lambda = wrapper_lambda*B**2.0_dp
    end function wrapper_lambda

    subroutine plot_lambda(theta_0, lambda_integral_analytic, lambda_integral)
        use myplot_module, only: myplot
        real(dp), dimension(:), intent(in) :: theta_0
        real(dp), dimension(:), intent(in) :: lambda_integral_analytic
        real(dp), dimension(:), intent(in) :: lambda_integral

        type(myplot) :: plt

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$\int \lambda d\varphi$ [1]", &
                            legend=.true.)
        call plt%add_plot(theta_0, &
                          lambda_integral_analytic, &
                          label="analytic", &
                          linestyle="r-")
        call plt%add_plot(theta_0, &
                          lambda_integral, &
                          label="numeric", &
                          linestyle="b--")
        call plt%show()
    end subroutine plot_lambda

end program test_lambda_helical_anti_sigma
