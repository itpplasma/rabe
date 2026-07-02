program test_shaing_callen_remainder_scaling
    use mock_field, only: mock_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use constants, only: dp, pi
    use utils, only: linspace, not_same

    use shaing_callen_remainder, only: get_non_omnigenous_remainder_pitch
    use shaing_callen_mod, only: calc_trapped_fraction
    use shaing_callen_remainder, only: calc_trapped_fraction_prime

    use plot_quantities, only: plot_non_omnigenous_remainder

    implicit none

    integer, parameter :: n_fieldlines = 50

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(flock_of_fieldlines_t) :: flock

    real(dp), parameter :: nfp = 4.0_dp, iota = 1.24_dp
    real(dp), parameter :: M_pol = -1.0_dp, N_tor = nfp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.15_dp
    type(mock_field_t) :: qs_field

    real(dp), parameter :: M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    real(dp) :: eps_1
    type(mock_perturbed_field_t) :: field

    integer, parameter :: n_eta = 400
    integer :: this
    integer, parameter :: n_eps = 8
    real(dp), dimension(n_eps) :: eps, remainder_pitch, remainder
    real(dp), parameter :: estimated_scaling_pitch = abs(1.42*M_pol/nfp)* &
                           (1.0_dp - 2.0_dp*abs(eps_0))
    real(dp), parameter :: estimated_scaling = (1.0_dp - 2.0_dp*abs(eps_0)) &
                           /sqrt(abs(eps_0))
    real(dp) :: abstol

    logical, parameter :: should_plot = .false.
    logical :: test_failed

    test_failed = .false.

    call qs_field%mock_field_init(M_pol, N_tor, B_0, eps_0)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call linspace(2.0_dp, real(n_eps, kind=dp), n_eps, eps)
    eps = 10.0_dp**(-eps)

    do this = 1, n_eps
        eps_1 = eps(this)*(B_0 + abs(eps_0))
        call field%mock_perturbed_field_init(qs_field, M_pol_pert, N_tor_pert, eps_1)
        call make_flock_from_labels(flock, &
                                    xi_0, &
                                    iota, &
                                    field, &
                                    M_pol, &
                                    N_tor, &
                                    nfp)

        remainder_pitch(this) = get_non_omnigenous_remainder_pitch(flock, &
                                                                   field, &
                                                                   n_eta)
        abstol = estimated_scaling_pitch*sqrt(abs(eps(this)))
        if (not_same(remainder_pitch(this), 0.0_dp, &
                     reltol_in=0.0_dp, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_get_non_omnigenous_remainder_scaling failed: ", &
                "pitch boundary term"
            print *, "omnigenity violation = ", eps(this)
            print *, "abs remainder = ", abs(remainder_pitch(this))
            print *, "expected remainder = ", abstol
            test_failed = .true.
        end if

        remainder(this) = calc_trapped_fraction(flock, field, n_eta)* &
                          M_pol/(M_pol*iota - N_tor)
        remainder(this) = 1.0_dp - &
                          calc_trapped_fraction_prime(flock, field, n_eta)/ &
                          remainder(this)

        ! We increase tolerance of ~20% as for the analytic estimate only the
        ! contribution of the above pitch boundary term remainder was considered
        ! and other linear terms where ignored
        abstol = estimated_scaling*sqrt(abs(eps(this)))*1.2_dp
        if (not_same(remainder(this), 0.0_dp, reltol_in=0.0_dp, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_get_non_omnigenous_remainder_scaling failed: ", &
                "1 - (M\iota - N)f_t'/(M f_t)"
            print *, "omnigenity violation = ", eps(this)
            print *, "abs remainder = ", abs(remainder(this))
            print *, "expected remainder = ", abstol
            test_failed = .true.
        end if
    end do

    if (should_plot) then
        call plot_non_omnigenous_remainder(eps, &
                                           remainder_pitch, &
                                           "pitch boundary term", &
                                           estimated_scaling_pitch)
        call plot_non_omnigenous_remainder(eps, &
                                           remainder, &
                                         "$1 - \frac{M\iota - N}{M}\frac{f_t'}{f_t}$", &
                                           estimated_scaling)
    end if

    if (test_failed) error stop

end program test_shaing_callen_remainder_scaling
