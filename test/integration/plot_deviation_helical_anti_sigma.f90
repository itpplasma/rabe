program plot_deviation_helical_anti_sigma
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_delta_eta
    use plot_quantities, only: plot_delta_A
    use plot_quantities, only: plot_deviation

    implicit none

    real(dp), parameter :: M_pol = 2.0_dp, N_tor = 10.0_dp
    character(len=*), parameter :: bc_filename = "input/helical_anti.bc"
    real(dp), parameter :: psi_edge = abs(-0.28274_dp)/(2.0_dp*pi) !Tm^2
    real(dp), parameter :: dr_dpsi = 1.0_dp/0.0661777_dp/psi_edge/100.0_dp
    real(dp), parameter :: J_pol_over_N_tor = -4.0_dp*1e6
    real(dp), parameter :: I_tor = 0.0_dp
    real(dp), parameter :: stor = 0.9999_dp, R = 8.0_dp
    real(dp), parameter :: eps_0 = -0.05, eps_1 = -0.0375
    real(dp), parameter :: eps_ratio = eps_1/abs(eps_0)
    real(dp), parameter :: delta_A_1 = 0.25_dp*eps_ratio*(1.0_dp + 6.0_dp*abs(eps_0))
    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 7e-7
    integer, parameter :: n_fieldlines = 20

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

    call plot_fieldlines_over_field(fieldlines, field, N_tor)
    call plot_delta_eta(fieldlines)
    call plot_delta_A(fieldlines, delta_A_1)

    call calc_deviation(fieldlines, field, deviation_A, deviation_B)

    covariant_factor = -2.0_dp*1e-7*(J_pol_over_N_tor*abs(N_tor) + I_tor*iota)
    off_factor_A = deviation_A*dr_dpsi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dpsi

    call plot_deviation(off_factor_A, off_factor_B)

end program plot_deviation_helical_anti_sigma
