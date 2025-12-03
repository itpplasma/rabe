program plot_helical_omnigenous
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_phi_max_over_xi_0
    use plot_quantities, only: plot_delta_A, plot_delta_eta
    use plot_quantities, only: plot_deviation_spectrum

    implicit none

    real(dp), parameter :: M_pol = 1.0_dp, N_tor = 4.0_dp
    character(len=*), parameter :: filename = "input/hsx_omnigenous.bc"
    real(dp), parameter :: stor = 0.25
    real(dp), parameter :: ds_dr = 0.190236_dp, sign_sqrtg = -1.0_dp
    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 8e-6
    integer, parameter :: n_fieldlines = 101

    real(dp) :: iota, nfp
    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

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

    call plot_fieldlines_over_field(fieldlines, field)
    call plot_phi_max_over_xi_0(fieldlines)
    call plot_delta_A(fieldlines)
    call plot_delta_eta(fieldlines)
    call plot_deviation_spectrum(fieldlines)

end program plot_helical_omnigenous
