program plot_deviation_goodman_squid
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_labels, only: get_labels
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation
    use fit_functions, only: S_A, S_B

    use plot_quantities, only: plot_deviation_spectrum
    use plot_quantities, only: plot_distribution_function
    use readers, only: read_column
    use plot_quantities, only: external_data_t
    use plot_quantities, only: plot_delta_eta_modes
    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_maxima_over_label
    use plot_quantities, only: plot_delta_eta
    use plot_quantities, only: plot_delta_A
    use plot_quantities, only: plot_deviation, external_data_t
    use plot_quantities, only: plot_phi_max_over_xi_0

    implicit none

    character(len=*), parameter :: bc_filename = "input/nautilus.bc"

    !------------------Taken from .bc file-------------------------------------!
    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 4.0_dp
    real(dp), parameter :: sign_sqrtg = -1.0_dp ! theta goes counter-clockwise
    !--------------------------------------------------------------------------!

    real(dp), parameter :: stor = 0.25_dp

    !------------------Taken from NEO-2 output---------------------------------!
    real(dp), parameter :: ds_dr = 0.00752119_dp*100.0_dp ! [1/m] called "avnabpsi"
    !--------------------------------------------------------------------------!

    real(dp) :: R
    real(dp) :: dr_dAphi

    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 8e-7
    integer, parameter :: max_n_fieldlines = 151

    real(dp), dimension(:), allocatable :: xi_0
    real(dp) :: iota, nfp
    real(dp) :: approx_iota
    integer :: n_fieldlines
    type(fieldline_t), dimension(:), allocatable :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    logical, parameter :: should_plot_others = .true.

    type(external_data_t) :: g_neo2
    real(dp), parameter :: nu_star = 6e-5
    real(dp), parameter :: scaling = 100.0_dp/5.884 !100/bmod0 -> g_(NEO-2) in [cm]

    call field%neo_field_init(bc_filename, stor)
    iota = field%iota
    nfp = field%nfp
    call get_labels(max_n_fieldlines, iota, M_pol, N_tor, nfp, xi_0, approx_iota)
    n_fieldlines = size(xi_0)
    allocate (fieldlines(n_fieldlines))

    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  approx_iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    if (should_plot_others) then
        g_neo2%label = "NEO-2: $\hat{g}_0 \frac{\mathrm{bmod0}}{100}$"
        call read_column("input/theta.dat", g_neo2%x, 1, 1)
        call read_column("input/gvpar0.dat", &
                         g_neo2%y, &
                         size(g_neo2%x, dim=1), &
                         size(g_neo2%x, dim=1)/8)
        g_neo2%y = g_neo2%y/scaling
        call plot_distribution_function(fieldlines, field, nu_star, g_neo2)
        call plot_deviation_spectrum(fieldlines)
        call plot_delta_eta_modes(fieldlines)
        call plot_fieldlines_over_field(fieldlines, field)
        call plot_maxima_over_label(fieldlines)
        call plot_delta_eta(fieldlines)
        call plot_delta_A(fieldlines)
        call plot_phi_max_over_xi_0(fieldlines)
    end if

    call calc_deviation(fieldlines, deviation_A, deviation_B)

    covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*approx_iota)
    dr_dAphi = 1.0_dp/(ds_dr*field%psi_tor_edge)*sign_sqrtg
    R = field%R
    off_factor_A = deviation_A*dr_dAphi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dAphi

    call plot_deviation(off_factor_A, &
                        off_factor_B)

end program plot_deviation_goodman_squid
