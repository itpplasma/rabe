program test_shaing_callen_against_quasi_symmetric
    use mock_field, only: mock_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use constants, only: dp, pi
    use utils, only: linspace
    use test_shaing_callen_mod, only: calc_quasi_symmetric_trapped_fraction
    use shaing_callen_mod, only: get_eta_integration_grid
    use shaing_callen_mod, only: calc_eta_integrand, calc_shaing_callen
    use shaing_callen_mod, only: eta_integrand_t

    implicit none

    real(dp) :: phi_tol = 7e-7
    integer, parameter :: n_fieldlines = 40

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp), parameter :: nfp = 10.0_dp, iota = 0.47_dp
    real(dp), parameter :: M_pol = 2.0_dp, N_tor = nfp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.01_dp
    type(mock_field_t) :: field

    integer :: this
    real(dp) :: trapped_fraction
    integer, parameter :: n_eta = 100
    real(dp), dimension(n_eta) :: eta_grid
    type(eta_integrand_t), dimension(n_fieldlines) :: eta_integrands
    real(dp) :: shaing_callen

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
    trapped_fraction = calc_quasi_symmetric_trapped_fraction(field, &
                                                             fieldlines(1)%eta_b)

    eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
    do this = 1, n_fieldlines
        call calc_eta_integrand(field, &
                                fieldlines(this), &
                                eta_grid, &
                                eta_integrands(this))
    end do
    shaing_callen = calc_shaing_callen(fieldlines, eta_integrands)

    if (test_failed) error stop

end program test_shaing_callen_against_quasi_symmetric
