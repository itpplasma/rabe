program test_shaing_callen_against_quasi_symmetric
    use mock_field, only: mock_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use test_shaing_callen_mod, only: calc_quasi_symmetric_trapped_fraction
    use shaing_callen_integration, only: get_eta_integration_grid
    use shaing_callen_mod, only: calc_shaing_callen
    use shaing_callen_mod, only: calc_shaing_callen_prime
    use shaing_callen_mod, only: eta_integrand_t, shaing_callen_t
    use shaing_callen_mod, only: calc_F, calc_alternative_F
    use shaing_callen_mod, only: calc_eta_integrand
    use test_calc_eta_integrand_mod, only: plot_eta_integrands

    implicit none

    real(dp), parameter :: reltol = 1e-6, abstol = 0.0_dp
    real(dp), parameter :: phi_tol = 7e-7
    integer, parameter :: n_fieldlines = 2

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp), parameter :: nfp = 10.0_dp, iota = 0.47_dp
    real(dp), parameter :: M_pol = 2.0_dp, N_tor = nfp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.01_dp
    type(mock_field_t) :: field

    integer :: this
    real(dp) :: qs_trapped_fraction
    real(dp) :: qs_modified_trapped_fraction
    integer, parameter :: n_eta = 10
    real(dp), dimension(n_eta) :: qs_F, F, altervative_F
    real(dp), dimension(n_eta) :: eta_grid
    type(shaing_callen_t) :: shaing_callen
    type(shaing_callen_t) :: shaing_callen_prime
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
    qs_trapped_fraction = calc_quasi_symmetric_trapped_fraction(field, &
                                                                fieldlines(1)%eta_b, &
                                                                n_eta)

    eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
    shaing_callen = calc_shaing_callen(field, fieldlines, eta_grid)
    shaing_callen_prime = calc_shaing_callen_prime(field, fieldlines, eta_grid)

    do this = 1, n_fieldlines
        call calc_eta_integrand(field, fieldlines(this), eta_grid, eta_integrands(this))
    end do
    call plot_eta_integrands(eta_integrands, fieldlines)

    if (not_same(qs_trapped_fraction, &
                 shaing_callen%trapped_fraction, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_shaing_callen_against_quasi_symmetric failed: ", &
            "qs_trapped_particle_fraction =/= trapped_particle_fraction"
        print *, "qs_trapped_particle_fraction = ", qs_trapped_fraction
        print *, "trapped_particle_fraction = ", shaing_callen%trapped_fraction
        print *, "relative error = ", &
            abs(1.0_dp - shaing_callen%trapped_fraction/qs_trapped_fraction)
        test_failed = .true.
    end if

    qs_modified_trapped_fraction = qs_trapped_fraction*M_pol/ &
                                   (M_pol*iota - N_tor)

    altervative_F = calc_alternative_F(field, fieldlines, eta_grid)
    F = calc_F(field, fieldlines, eta_grid)
    qs_F = M_pol/(M_pol*iota - N_tor)*eta_grid

    if (not_same(qs_F, &
                 F, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_shaing_callen_against_quasi_symmetric failed: ", &
            "qs_F(eta) =/= F(eta)"
        print *, "qs_F(eta) = ", qs_F
        print *, "F(eta) = ", F
        print *, "relative error = ", abs(1.0_dp - F/qs_F)
        test_failed = .true.
    end if

    if (not_same(qs_F, &
                 altervative_F, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_shaing_callen_against_quasi_symmetric failed: ", &
            "qs_F(eta) =/= altervative_F(eta)"
        print *, "qs_F(eta) = ", qs_F
        print *, "altervative_F(eta) = ", F
        print *, "relative error = ", abs(1.0_dp - altervative_F/qs_F)
        test_failed = .true.
    end if

    if (not_same(qs_modified_trapped_fraction, &
                 shaing_callen%modified_trapped_fraction, &
                 reltol_in=reltol, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_shaing_callen_against_quasi_symmetric failed: ", &
            "qs_modified_trapped_fraction =/= modified_trapped_fraction"
        print *, "qs_modified_trapped_fraction = ", qs_modified_trapped_fraction
        print *, "trapped_particle_fraction = ", shaing_callen%modified_trapped_fraction
        print *, "relative error = ", &
            abs(1.0_dp - shaing_callen%modified_trapped_fraction/ &
                qs_modified_trapped_fraction)
        test_failed = .true.
    end if

    if (test_failed) error stop

end program test_shaing_callen_against_quasi_symmetric
