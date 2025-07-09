program plot_deviation_goodman_squid
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

    real(dp), parameter :: phi_tol = 5e-7
    integer, parameter :: n_fieldlines = 50

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    logical, parameter :: should_plot_others = .true.

    call field%neo_field_init(bc_filename, stor)
    iota = field%iota
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
        call plot_maxima_over_label(fieldlines)
        call plot_delta_eta(fieldlines)
        call plot_delta_A(fieldlines)
    end if

    call calc_deviation(fieldlines, field, deviation_A, deviation_B)

    covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*iota)
    dr_dAphi = 1.0_dp/(ds_dr*field%psi_tor_edge)*sign_sqrtg
    R = field%R
    off_factor_A = deviation_A*dr_dAphi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dAphi

    call plot_deviation(off_factor_A, &
                        off_factor_B)

end program plot_deviation_goodman_squid
