program plot_deviation_landreman_paul_helical
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation
    use misc, only: S_A, S_B

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_maxima_over_label
    use plot_quantities, only: plot_delta_eta
    use plot_quantities, only: plot_delta_A
    use plot_quantities, only: plot_deviation, external_data_t

    implicit none

    character(len=*), parameter :: bc_filename = "input/landreman_paul_helical.bc"

    !------------------Taken from .bc file-------------------------------------!
    real(dp), parameter :: M_pol = 1.0_dp, N_tor = -4.0_dp
    real(dp), parameter :: sign_sqrtg = -1.0_dp ! theta goes counter-clockwise
    !--------------------------------------------------------------------------!

    real(dp), parameter :: stor = 0.60_dp

    !------------------Taken from NEO-2 output---------------------------------!
    real(dp), parameter :: ds_dr = 0.00852345_dp*100.0_dp ! [1/m] called "avnabpsi"
    integer, parameter :: n_neo2 = 18
    real(dp), dimension(n_neo2), parameter :: nu_star_neo2 = (/3e-09, &
                                                               1e-08, &
                                                               3e-08, &
                                                               1e-07, &
                                                               3e-07, &
                                                               4e-07, &
                                                               7.1e-07, &
                                                               1e-06, &
                                                               1.2e-06, &
                                                               2e-06, &
                                                               3e-06, &
                                                               6e-06, &
                                                               1e-05, &
                                                               1.1e-05, &
                                                               1.9e-05, &
                                                               3e-05, &
                                                               5.8e-05, &
                                                               0.0001/)
    real(dp), dimension(n_neo2), parameter :: lambda_neo2 = (/7.29201_dp, &
                                                              6.15848_dp, &
                                                              3.69955_dp, &
                                                              1.70273_dp, &
                                                              -0.293197_dp, &
                                                              -0.814098_dp, &
                                                              -1.69767_dp, &
                                                              -2.10485_dp, &
                                                              -2.28529_dp, &
                                                              -2.67295_dp, &
                                                              -2.8802_dp, &
                                                              -3.09469_dp, &
                                                              -3.18068_dp, &
                                                              -3.19218_dp, &
                                                              -3.23892_dp, &
                                                              -3.25972_dp, &
                                                              -3.2713_dp, &
                                                              -3.27004_dp/)
    real(dp), parameter :: shaing_callen_limit = -3.3268_dp
    type(external_data_t) :: lambda_off_neo2
    !--------------------------------------------------------------------------!

    real(dp) :: R
    real(dp) :: dr_dAphi

    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 5e-7
    integer, parameter :: n_fieldlines = 50

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota, nfp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    logical, parameter :: should_plot_others = .true.

    call field%neo_field_init(bc_filename, stor)
    iota = field%iota
    nfp = field%nfp
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, &
                                  theta_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    if (should_plot_others) then
        call plot_fieldlines_over_field(fieldlines, field, N_tor)
        call plot_maxima_over_label(fieldlines)
        call plot_delta_eta(fieldlines)
        call plot_delta_A(fieldlines)
    end if

    call calc_deviation(fieldlines, deviation_A, deviation_B)

    covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*iota)
    dr_dAphi = 1.0_dp/(ds_dr*field%psi_tor_edge)*sign_sqrtg
    R = field%R
    off_factor_A = deviation_A*dr_dAphi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dAphi

    allocate (lambda_off_neo2%x(n_neo2), lambda_off_neo2%y(n_neo2))
    lambda_off_neo2%label = "NEO-2: $\lambda_{bB} - \lambda^\mathrm{SC}$"
    lambda_off_neo2%x = nu_star_neo2
    lambda_off_neo2%y = lambda_neo2 - shaing_callen_limit
    lambda_off_neo2%y = lambda_off_neo2%y + 3e-2
    call plot_deviation(off_factor_A, &
                        off_factor_B, &
                        lambda_off_external=lambda_off_neo2)
    deallocate (lambda_off_neo2%x, lambda_off_neo2%y)

end program plot_deviation_landreman_paul_helical
