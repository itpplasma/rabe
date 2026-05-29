program test_set_fieldline_labels
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use field_checks, only: suspect_omnigenous_origin_not_minimum
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_fieldlines_over_field_chi_xi
    use plot_quantities, only: plot_delta_A
    use plot_quantities, only: plot_maxima_over_label

    implicit none

    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.01, eps_1 = -0.0001_dp
    type(anti_sigma_field_t) :: base_field
    real(dp), parameter :: B_pert = 0.0006_dp, N_pert = 2.0_dp, M_pert = 1.0_dp
    type(mock_perturbed_field_t) :: field

    real(dp), parameter :: chi_tol = 1.2e-4
    integer, parameter :: n_fieldlines = 100
    real(dp), dimension(2) :: phi_max

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp

    type(flock_of_fieldlines_t) :: flock

    integer, parameter :: n_cases = 4
    real(dp), dimension(n_cases) :: M_pols, N_tors, iotas, nfps
    real(dp) :: M_pol, N_tor, iota, nfp

    real(dp) :: chi_0
    real(dp) :: expected_chi_0
    real(dp) :: dummy_A, dummy_B
    integer :: case, current
    logical :: test_failed

    logical, parameter :: should_plot = .false.

    M_pols = [1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp]
    N_tors = [1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp]
    iotas = [-2.3_dp, 1.0_dp, -1.5_dp, -0.8_dp]
    nfps = [1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp]

    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)
    expected_chi_0 = 0.0_dp

    test_failed = .false.

    do case = 1, n_cases
        M_pol = M_pols(case)
        N_tor = N_tors(case)
        iota = iotas(case)
        nfp = nfps(case)
        print *, "--------------------------------------------------------------"
        print *, "case: M_pol=", M_pol, " N_tor=", N_tor, "iota=", iota

        call base_field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
        call field%mock_perturbed_field_init(base_field, M_pert, N_pert, B_pert)
        if (suspect_omnigenous_origin_not_minimum(field, M_pol, N_tor)) then
            print *, "error: The origin of the IDEAL omnigenous configuration"
            print *, "(theta=phi=0) must be a global and local minimum!"
            print *, "Origin of provided field suggests that this is not the case!"
            error stop
        end if
        do current = 1, n_fieldlines
            chi_0 = M_pol*(N_tor*xi_0(current)/nfp) - N_tor*(M_pol*xi_0(current)/nfp)
            if (not_same(chi_0, expected_chi_0, abstol_in=chi_tol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_set_fieldline_labels_along_chi_min failed: chi_0"
                print *, "found: ", chi_0
                print *, "expected: ", expected_chi_0
                test_failed = .true.
            end if
        end do
        call make_flock_of_fieldlines(flock, xi_0, iota, field, M_pol, N_tor, &
                                      nfp)
        if (should_plot) then
            call plot_fieldlines_over_field_chi_xi(flock%fieldlines, field, &
                                                   flock%M_pol, flock%N_tor, flock%nfp)
            call plot_maxima_over_label(flock%fieldlines, flock%iota_p)
            call plot_fieldlines_over_field(flock%fieldlines, field)
            call plot_delta_A(flock%fieldlines)
        end if
        call calc_deviation(flock, dummy_A, dummy_B)
    end do

    if (test_failed) error stop

end program test_set_fieldline_labels
