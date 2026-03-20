program plot_deviation_drevlak_axi
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_labels, only: get_theta_0
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation
    use fit_functions, only: S_A, S_B

    use plot_quantities, only: plot_deviation_spectrum
    use plot_quantities, only: plot_delta_eta_modes
    use plot_quantities, only: plot_spectra
    use plot_quantities, only: plot_B_along_fieldline
    use plot_quantities, only: plot_fieldlines_over_field_chi_xi
    use plot_quantities, only: plot_maxima_over_label
    use plot_quantities, only: plot_delta_eta
    use plot_quantities, only: plot_delta_A
    use plot_quantities, only: plot_deviation, external_data_t

    implicit none

    character(len=*), parameter :: bc_filename = "input/landreman_paul_qa.bc"

    !------------------Taken from .bc file-------------------------------------!
    real(dp), parameter :: M_pol = 1.0_dp, N_tor = 0.0_dp
    real(dp), parameter :: sign_sqrtg = -1.0_dp ! theta goes counter-clockwise
    !--------------------------------------------------------------------------!

    real(dp), parameter :: stor = 0.496_dp

    !------------------Taken from NEO-2 output---------------------------------!
    real(dp), parameter :: ds_dr = 1.0_dp ! [1/m] called "avnabpsi"
    !--------------------------------------------------------------------------!

    real(dp) :: R
    real(dp) :: dr_dAphi

    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 1e-5
    integer, parameter :: max_n_fieldlines = 400

    real(dp), dimension(:), allocatable :: theta_0
    real(dp) :: iota, nfp
    real(dp) :: approx_iota
    integer :: n_fieldlines
    type(fieldline_t), dimension(:), allocatable :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    logical, parameter :: should_plot_others = .true.

    call field%neo_field_init(bc_filename, stor)
    iota = field%iota
    nfp = field%nfp
    call get_theta_0(max_n_fieldlines, iota, M_pol, N_tor, nfp, theta_0, approx_iota)
    n_fieldlines = size(theta_0)
    allocate (fieldlines(n_fieldlines))

    call make_flock_of_fieldlines(fieldlines, &
                                  theta_0, &
                                  approx_iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    if (should_plot_others) then
        call plot_delta_eta(fieldlines)
        call plot_delta_A(fieldlines)
        call plot_spectra(fieldlines)
        call plot_fieldlines_over_field_chi_xi(fieldlines, field, M_pol, N_tor, nfp)
        call plot_deviation_spectrum(fieldlines)
        call plot_maxima_over_label(fieldlines)
    end if

    call calc_deviation(fieldlines, deviation_A, deviation_B)

    covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*iota)
    dr_dAphi = 1.0_dp/(ds_dr*field%psi_tor_edge)*sign_sqrtg
    R = field%R
    off_factor_A = deviation_A*dr_dAphi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dAphi

    call plot_deviation(off_factor_A, &
                        off_factor_B)

end program plot_deviation_drevlak_axi
