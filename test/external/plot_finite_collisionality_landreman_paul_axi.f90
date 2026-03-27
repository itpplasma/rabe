program plot_finite_collisionality_landreman_paul_axi
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_labels, only: get_labels
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation
    use surface_average_mod, only: surface_average_t, calc_surface_averages
    use coefficients, only: calc_finite_boundary_layer_correction
    use shaing_callen_mod, only: calc_trapped_fraction

    use plot_quantities, only: plot_asymptotic_model

    implicit none

    character(len=*), parameter :: bc_filename = "input/landreman_paul_qa.bc"

    !------------------Taken from .bc file-------------------------------------!
    real(dp), parameter :: M_pol = 1.0_dp, N_tor = 0.0_dp
    real(dp), parameter :: sign_sqrtg = -1.0_dp ! theta goes counter-clockwise
    !--------------------------------------------------------------------------!

    real(dp), parameter :: stor = 0.496_dp

    real(dp) :: R
    real(dp) :: dr_dAtheta

    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 1e-5
    integer, parameter :: max_n_fieldlines = 400

    real(dp), dimension(:), allocatable :: xi_0
    real(dp) :: iota, nfp
    real(dp) :: approx_iota
    integer :: n_fieldlines
    type(fieldline_t), dimension(:), allocatable :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    integer, parameter :: n_eta = 100
    real(dp) :: trapped_fraction
    real(dp) :: lambda_SC
    real(dp) :: helical_factor
    real(dp) :: I_ref_hat
    real(dp) :: finite_col_correction
    type(surface_average_t) :: average

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

    call calc_deviation(fieldlines, deviation_A, deviation_B)

    covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*iota)
    call calc_surface_averages(fieldlines, average)
    dr_dAtheta = sign_sqrtg*sign(1.0_dp, field%psi_tor_edge)/average%sqrt_g11
    R = field%R
    off_factor_A = deviation_A*dr_dAtheta*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5_dp*R*pi*dr_dAtheta

    trapped_fraction = calc_trapped_fraction(field, fieldlines, n_eta)
    helical_factor = (field%B_phi_covariant*M_pol + field%B_theta_covariant*N_tor)/ &
                     (M_pol*iota - N_tor)
    lambda_SC = helical_factor*dr_dAtheta*trapped_fraction

    finite_col_correction = calc_finite_boundary_layer_correction(fieldlines, &
                                                                  R, &
                                                                  dr_dAtheta, &
                                                              field%B_theta_covariant, &
                                                                  field%B_phi_covariant)

    call plot_asymptotic_model(off_factor_A, &
                               off_factor_B, &
                               lambda_SC, &
                               finite_col_correction)

end program plot_finite_collisionality_landreman_paul_axi
