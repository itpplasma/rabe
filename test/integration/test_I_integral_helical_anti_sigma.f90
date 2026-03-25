program test_I_integral_helical_anti_sigma
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use fieldline_labels, only: get_labels
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_delta_eta
    use plot_quantities, only: plot_I
    use plot_quantities, only: plot_I_factor
    use plot_quantities, only: plot_delta_A

    implicit none

    real(dp), parameter :: M_pol = 2.0_dp, N_tor = 10.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.05, eps_1 = -0.00375_dp
    real(dp), parameter :: factor = sqrt(abs(eps_0)/(1.0_dp + abs(eps_0))) &
                           /B_0**2.0_dp/N_tor*sqrt(32.0_dp)
    real(dp), parameter :: I_0_analytic = factor*(1.0_dp + abs(eps_0)*2.0_dp/3.0_dp)
    real(dp), parameter :: eps_ratio = eps_1/abs(eps_0)
    real(dp), parameter :: I_1_analytic = factor*(-0.5_dp*eps_ratio)* &
                           (1.0_dp + 6.0_dp*abs(eps_0))
    real(dp), parameter :: B_max = B_0*(1.0_dp + abs(eps_0))
    real(dp), parameter :: delta_A_1 = 0.25_dp*eps_ratio*(1.0_dp + 6.0_dp*abs(eps_0))
    real(dp), parameter :: psi_edge = abs(-0.28274_dp)/(2.0_dp*pi) !Tm^2
    real(dp), parameter :: dr_dpsi = 1.0_dp/0.0661777_dp/psi_edge/100.0_dp
    real(dp), parameter :: J_pol_over_N_tor = -4.0_dp*1e6
    real(dp), parameter :: I_tor = 0.0_dp
    real(dp), parameter :: R = 8.0_dp
    type(anti_sigma_field_t) :: field
    real(dp), dimension(2), parameter :: chi_max = [+pi, -pi]

    real(dp), parameter :: phi_tol = 1e-6
    real(dp), parameter :: B_max_reltol = phi_tol**2.0_dp*N_tor**2.0
    real(dp), parameter :: I_mean_reltol = eps_0*eps_0*3.0_dp
    real(dp), parameter :: I_amplitude_reltol = max(eps_ratio*eps_ratio*5.0_dp, &
                                                    eps_0*eps_0*10.0_dp)
    integer, parameter :: max_n_fieldlines = 20
    real(dp), dimension(2) :: phi_max

    real(dp), dimension(:), allocatable :: xi_0
    real(dp), parameter :: iota = 0.0_dp
    real(dp), parameter :: nfp = max(1.0_dp, N_tor)
    real(dp) :: approx_iota
    integer :: n_fieldlines
    type(fieldline_t), dimension(:), allocatable :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    integer :: current
    real(dp), dimension(:), allocatable :: I_shifted
    real(dp) :: I_mean, I_amplitude
    logical :: test_failed

    logical, parameter :: should_plot = .false.

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call get_labels(max_n_fieldlines, iota, M_pol, N_tor, nfp, xi_0, approx_iota)
    n_fieldlines = size(xi_0)
    allocate (fieldlines(n_fieldlines))
    allocate (I_shifted(n_fieldlines))
    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  approx_iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    test_failed = .false.
    do current = 1, n_fieldlines
        phi_max = (M_pol*fieldlines(current)%theta_0 - chi_max)/N_tor
        if (not_same(fieldlines(current)%phi_max(1), &
                     phi_max(1), &
                     abstol_in=2.0_dp*phi_tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_I_integral_helical_anti_sigma failed: left phi_max"
            print *, "found: ", fieldlines(current)%phi_max(1)
            print *, "analytic: ", phi_max(1)
            test_failed = .true.
        end if
        if (not_same(fieldlines(current)%phi_max(2), &
                     phi_max(2), &
                     abstol_in=2.0_dp*phi_tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_I_integral_helical_anti_sigma failed: right phi_max"
            print *, "found: ", fieldlines(current)%phi_max(2)
            print *, "analytic: ", phi_max(2)
            test_failed = .true.
        end if
     if (not_same(1.0_dp/fieldlines(current)%eta_b, B_max, reltol_in=B_max_reltol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_I_integral_helical_anti_sigma failed: B_max"
            print *, "found: ", 1.0_dp/fieldlines(current)%eta_b
            print *, "analytic: ", B_max
            print *, "relative error: ", fieldlines(current)%eta_b*B_max - 1.0_dp
            test_failed = .true.
        end if
    end do

    I_mean = sum(fieldlines%integral_lambda_b_over_B_squared)/n_fieldlines
    if (not_same(I_mean, I_0_analytic, reltol_in=I_mean_reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_I_integral_helical_anti_sigma failed: I_mean"
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
        print *, "test_I_integral_helical_anti_sigma failed: I_amplitude"
        print *, "found: ", I_amplitude
        print *, "analytic: ", I_1_analytic
        print *, "relative error: ", I_amplitude/I_1_analytic - 1.0_dp
        test_failed = .true.
    end if

    if (test_failed) error stop

    if (should_plot) then
        call plot_fieldlines_over_field(fieldlines, field)
        call plot_delta_eta(fieldlines)
        call plot_I_factor(fieldlines, eps_0, eps_1)
        call plot_I(fieldlines, I_0_analytic, I_1_analytic, eps_0, eps_1)
        call plot_delta_A(fieldlines, delta_A_1)
    end if

end program test_I_integral_helical_anti_sigma
