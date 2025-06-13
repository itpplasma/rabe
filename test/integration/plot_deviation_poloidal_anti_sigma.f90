program plot_deviation_poloidal_anti_sigma
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation
    use misc, only: S_A, S_B

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_delta_eta
    use plot_quantities, only: plot_delta_A
    use plot_quantities, only: plot_deviation

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    character(len=*), parameter :: bc_filename = "input/poloidal_anti.bc"
    real(dp), parameter :: psi_edge = abs(-0.00785398_dp)/(2.0_dp*pi) !Tm^2
    real(dp), parameter :: R = 1.00_dp
    real(dp), parameter :: J_pol_over_N_tor = -5.0_dp*1e6
    real(dp), parameter :: I_tor = 0.0_dp

    real(dp), parameter :: stor = 0.25781_dp
    real(dp), parameter :: ds_dr = 0.387524_dp*100.0_dp !1/m
    real(dp), parameter :: dr_dpsi = 1.0_dp/(ds_dr*psi_edge)

    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = 0.125, eps_1 = 0.05_dp
    real(dp), parameter :: delta_B_1 = 0.0_dp*1e-4
    real(dp), parameter :: eps_ratio = eps_1/abs(eps_0)
    real(dp), parameter :: delta_A_1 = 0.25_dp*eps_ratio*(1.0_dp + 6.0_dp*abs(eps_0))
    real(dp), parameter :: B_max = B_0*(1.0_dp + eps_0) + abs(delta_B_1)
    real(dp), parameter :: delta_eta_1 = delta_B_1/B_max**2.0_dp
    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 4e-6
    integer, parameter :: n_fieldlines = 20

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.1618_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp) :: deviation_A, deviation_B
    integer, parameter :: n_points = 100
    real(dp), dimension(n_points) :: nu_star
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    real(dp), parameter :: iota_p = iota*pi/(N_tor - iota*M_pol)
    real(dp) :: off_factor_A_analytic, off_factor_B_analytic

    logical, parameter :: should_plot_others = .false.

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

    if (should_plot_others) then
        call plot_fieldlines_over_field(fieldlines, field, N_tor)
        call plot_delta_eta(fieldlines, delta_eta_1, fieldlines(1)%iota_p)
        call plot_delta_A(fieldlines, delta_A_1)
    end if

    call calc_deviation(fieldlines, field, deviation_A, deviation_B)

    covariant_factor = -2.0_dp*1e-7*(J_pol_over_N_tor*abs(N_tor) + I_tor*iota)
    off_factor_A = deviation_A*dr_dpsi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dpsi

    off_factor_A_analytic = -dr_dpsi*sqrt(2.0_dp*pi)*R*B_0*abs(eps_1)**2.0_dp &
                            *N_tor*S_A(iota_p)/(2.0_dp*abs(eps_0))**0.75_dp &
                            /(N_tor - iota*M_pol)**1.5_dp*0.5_dp &
                            *(1.0_dp + 6.0_dp*abs(eps_0))

    off_factor_B_analytic = eps_1*dr_dpsi*R*pi*delta_B_1*S_B(iota_p)*0.25_dp

    call plot_deviation(off_factor_A, &
                        off_factor_B, &
                        off_factor_A_analytic, &
                        off_factor_B_analytic)

end program plot_deviation_poloidal_anti_sigma
