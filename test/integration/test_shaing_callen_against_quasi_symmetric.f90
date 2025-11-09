program test_shaing_callen_against_quasi_symmetric
    use mock_field, only: mock_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use test_shaing_callen_mod, only: test_F_against_quasi_symmetric
    use test_shaing_callen_mod, only: test_alternative_F_against_quasi_symmetric
    use test_shaing_callen_mod, only: test_trapped_fractions_against_quasi_symmetric
    use test_shaing_callen_mod, only: test_calc_avg_normalized_B_squared_dphimax_dxi0
    use test_shaing_callen_mod, only: test_calc_avg_normalized_lambda_dphimax_dxi0
    use test_shaing_callen_mod, only: test_get_non_omnigenous_remainder

    use shaing_callen_integration, only: get_eta_integration_grid
    use shaing_callen_mod, only: eta_integrand_t
    use shaing_callen_mod, only: calc_eta_integrand
    use test_calc_eta_integrand_mod, only: plot_eta_integrands
    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_phi_max_over_xi_0

    implicit none

    real(dp), parameter :: phi_tol = 7e-7
    integer, parameter :: n_fieldlines = 6

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp), parameter :: nfp = 10.0_dp, iota = 0.47_dp
    real(dp), parameter :: M_pol = 2.0_dp, N_tor = nfp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.01_dp
    type(mock_field_t) :: field

    logical, parameter :: should_plot = .false.
    integer :: this
    integer, parameter :: n_eta = 10
    real(dp), dimension(n_eta) :: eta_grid
    type(eta_integrand_t), dimension(n_fieldlines) :: eta_integrands

    logical :: test_failed

    test_failed = .false.

    call field%mock_field_init(M_pol, N_tor, B_0, eps_0)

    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    if (should_plot) then
        eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
        do this = 1, n_fieldlines
            call calc_eta_integrand(field, &
                                    fieldlines(this), &
                                    eta_grid, &
                                    eta_integrands(this))
        end do
        call plot_eta_integrands(eta_integrands, fieldlines)
        call plot_fieldlines_over_field(fieldlines, field)
        call plot_phi_max_over_xi_0(fieldlines)
    end if

    call test_F_against_quasi_symmetric(field, fieldlines, test_failed)
    call test_alternative_F_against_quasi_symmetric(field, fieldlines, test_failed)
    call test_trapped_fractions_against_quasi_symmetric(field, fieldlines, test_failed)
    call test_calc_avg_normalized_B_squared_dphimax_dxi0(fieldlines, test_failed)
    call test_calc_avg_normalized_lambda_dphimax_dxi0(field, fieldlines, test_failed)
    call test_get_non_omnigenous_remainder(field, fieldlines, test_failed)

    if (test_failed) error stop

end program test_shaing_callen_against_quasi_symmetric
