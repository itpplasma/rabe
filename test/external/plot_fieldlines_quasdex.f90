program plot_fieldlines_quasdex
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_phi_max_over_xi_0
    use test_shaing_callen_mod, only: test_calc_avg_normalized_lambda_dphimax_dxi0
    use test_shaing_callen_mod, only: test_calc_avg_normalized_B_squared_dphimax_dxi0

    implicit none

    character(len=*), parameter :: bc_filename = "input/quasdex.bc"

    !------------------Taken from .bc file-------------------------------------!
    real(dp), parameter :: M_pol = 1.0_dp, N_tor = 0.0_dp
    !--------------------------------------------------------------------------!

    real(dp), parameter :: stor = 0.2503_dp

    type(neo_field_t) :: field

    integer, parameter :: n_fieldlines = 101

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota, nfp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    logical :: test_failed

    test_failed = .false.

    call field%neo_field_init(bc_filename, stor)
    iota = field%iota
    nfp = field%nfp
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp)

    call plot_fieldlines_over_field(fieldlines, field)
    call plot_phi_max_over_xi_0(fieldlines)

end program plot_fieldlines_quasdex
