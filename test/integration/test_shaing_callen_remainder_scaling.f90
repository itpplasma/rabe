program test_shaing_callen_remainder_scaling
    use mock_field, only: mock_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use constants, only: dp, pi
    use utils, only: linspace, not_same

    use shaing_callen_remainder, only: get_non_omnigenous_remainder_pitch

    use plot_quantities, only: plot_non_omnigenous_remainder

    implicit none

    real(dp), parameter :: phi_tol = 7e-7
    integer, parameter :: n_fieldlines = 50

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    real(dp), parameter :: nfp = 4.0_dp, iota = 1.24_dp
    real(dp), parameter :: M_pol = -1.0_dp, N_tor = nfp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.3_dp
    type(mock_field_t) :: qs_field

    real(dp), parameter :: M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    real(dp) :: eps_1
    type(mock_perturbed_field_t) :: field

    integer, parameter :: n_eta = 400
    integer :: this
    integer, parameter :: n_eps = 8
    real(dp), dimension(n_eps) :: eps, remainder
    real(dp), parameter :: abstol_factor = 2.46e-1
    real(dp) :: abstol

    logical, parameter :: should_plot = .false.
    logical :: test_failed

    test_failed = .false.

    call qs_field%mock_field_init(M_pol, N_tor, B_0, eps_0)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call linspace(1.0_dp, real(n_eps, kind=dp), n_eps, eps)
    eps = 10.0_dp**(-eps)

    do this = 1, n_eps
        eps_1 = eps(this)*(B_0 + abs(eps_0))
        call field%mock_perturbed_field_init(qs_field, M_pol_pert, N_tor_pert, eps_1)
        call make_flock_of_fieldlines(fieldlines, &
                                      xi_0, &
                                      iota, &
                                      field, &
                                      M_pol, &
                                      N_tor, &
                                      nfp, &
                                      phi_tol)

        remainder(this) = get_non_omnigenous_remainder_pitch(field, &
                                                             fieldlines, &
                                                             n_eta)
        abstol = abstol_factor*sqrt(abs(eps(this)))
        if (not_same(remainder(this), 0.0_dp, reltol_in=0.0_dp, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_get_non_omnigenous_remainder_scaling failed:"
            print *, "omnigenity violation = ", eps(this)
            print *, "abs remainder = ", abs(remainder(this))
            print *, "expected remainder = ", abstol
            test_failed = .true.
        end if
    end do

    if (should_plot) then
        call plot_non_omnigenous_remainder(eps, remainder)
    end if

    if (test_failed) error stop

end program test_shaing_callen_remainder_scaling
