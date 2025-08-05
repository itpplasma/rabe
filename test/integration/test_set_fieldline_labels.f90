program test_set_fieldline_labels
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: set_fieldline_labels_along_chi_min
    use make_fieldline, only: make_flock_of_fieldlines

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_fieldlines_over_field_chi_xi
    use plot_quantities, only: plot_delta_A

    implicit none

    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.01, eps_1 = -0.0001_dp
    type(anti_sigma_field_t) :: field

    real(dp), parameter :: phi_tol = 2e-6

    integer, parameter :: n_fieldlines = 20
    real(dp), dimension(2) :: phi_max

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp

    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    integer, parameter :: n_cases = 4
    real(dp), dimension(n_cases) :: M_pols, N_tors, iotas
    real(dp) :: M_pol, N_tor, iota

    real(dp) :: chi_0
    real(dp) :: expected_chi_0
    integer :: case, current
    logical :: test_failed

    logical, parameter :: should_plot = .false.

    M_pols = (/1.0_dp, 1.0_dp, 0.0_dp, -2.0_dp/)
    N_tors = (/1.0_dp, 0.0_dp, 2.0_dp, -1.0_dp/)
    iotas = (/2.5_dp, 100.0_dp, -1.5_dp, 0.0_dp/)

    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    test_failed = .false.

    do case = 1, n_cases
        M_pol = M_pols(case)
        N_tor = N_tors(case)
        iota = iotas(case)
        fieldlines%xi_0 = xi_0
        expected_chi_0 = 0.0_dp

        fieldlines%iota = iota
        call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
        call set_fieldline_labels_along_chi_min(field, M_pol, N_tor, fieldlines, &
                                                phi_tol)
        do current = 1, n_fieldlines
            chi_0 = M_pol*fieldlines(current)%theta_0 - N_tor*fieldlines(current)%phi_0
            if (not_same(chi_0, expected_chi_0, abstol_in=2.0_dp*phi_tol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_set_fieldline_labels_along_chi_min failed: chi_0"
                print *, "theta_0: ", fieldlines(current)%theta_0
                print *, "phi_0: ", fieldlines(current)%phi_0
                print *, "found: ", chi_0
                print *, "expected: ", expected_chi_0
                test_failed = .true.
            end if
        end do
        if (should_plot) then
            call make_flock_of_fieldlines(fieldlines, xi_0, iota, field, M_pol, N_tor, &
                                          phi_tol=phi_tol)
            call plot_fieldlines_over_field_chi_xi(fieldlines, field, M_pol, N_tor)
            call plot_fieldlines_over_field(fieldlines, field, N_tor)
            call plot_delta_A(fieldlines)
        end if
    end do

    if (test_failed) error stop

end program test_set_fieldline_labels
