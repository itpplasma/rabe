program plot_helical_omnigenous
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_phi_max_over_xi_0
    use plot_quantities, only: plot_delta_A, plot_delta_eta

    use deviation, only: calc_deviation

    use shaing_callen_mod, only: shaing_callen_t
    use shaing_callen_integration, only: get_eta_integration_grid
    use shaing_callen_mod, only: calc_shaing_callen
    use shaing_callen_mod, only: get_non_omnigenous_remainder

    implicit none

    real(dp), parameter :: M_pol = 1.0_dp, N_tor = 4.0_dp
    character(len=*), parameter :: filename = "input/hsx_omnigenous.bc"
    real(dp), parameter :: stor = 0.25
    real(dp), parameter :: ds_dr = 0.190236_dp, sign_sqrtg = -1.0_dp
    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 8e-6
    integer, parameter :: n_fieldlines = 101

    real(dp) :: iota, nfp, dr_dpsi, dr_dAtheta, covariant_factor
    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp) :: deviation_A, deviation_B, R, C_A, C_B

    integer, parameter :: n_eta = 100
    real(dp), dimension(n_eta) :: eta_grid
    type(shaing_callen_t) :: shaing_callen
    real(dp) :: lambda_SC, remainder

    logical, parameter :: should_plot = .false.
    logical, parameter :: calc_bootstrap_and_stop = .false.

    call field%neo_field_init(filename, stor)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)
    iota = field%iota
    nfp = field%nfp

    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    if (should_plot) then
        call plot_fieldlines_over_field(fieldlines, field)
        call plot_phi_max_over_xi_0(fieldlines)
        call plot_delta_A(fieldlines)
        call plot_delta_eta(fieldlines)
    end if

    if (calc_bootstrap_and_stop) then
        call calc_deviation(fieldlines, deviation_A, deviation_B)

        covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*iota)
        dr_dpsi = 1.0_dp/(ds_dr*100.0_dp)/field%psi_tor_edge
        dr_dAtheta = dr_dpsi*sign_sqrtg
        R = field%R
        C_A = deviation_A*dr_dAtheta*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
        C_B = deviation_B*0.5*R*pi*dr_dAtheta
        print *, "1/sqrt(nu_star) factor: ", C_A
        print *, "1/nu_star factor: ", C_B

        eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
        shaing_callen = calc_shaing_callen(field, fieldlines, eta_grid)
        lambda_SC = shaing_callen%modified_trapped_fraction*covariant_factor - &
                    shaing_callen%trapped_fraction*field%B_theta_covariant
        lambda_SC = lambda_SC*dr_dAtheta
        remainder = get_non_omnigenous_remainder(field, fieldlines, eta_grid)
        remainder = remainder*covariant_factor*dr_dAtheta
        print *, "lambda_SC: ", lambda_SC
        print *, "non-omnigneous remainder: ", remainder
        error stop
    end if

end program plot_helical_omnigenous
