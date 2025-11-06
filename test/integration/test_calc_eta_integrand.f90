program test_calc_eta_integrand
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use constants, only: dp, pi
    use utils, only: linspace
    use shaing_callen_mod, only: calc_eta_integrand, eta_integrand_t
    use shaing_callen_integration, only: get_eta_integration_grid
    use test_calc_eta_integrand_mod, only: average_eta_integrands_sum_equal_F
    use test_calc_eta_integrand_mod, only: plot_eta_integrands

    implicit none

    real(dp) :: phi_tol = 7e-7
    integer, parameter :: n_fieldlines = 40

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines
    type(eta_integrand_t), dimension(n_fieldlines) :: eta_integrands

    real(dp), parameter :: nfp = 10.0_dp, iota = 0.47_dp
    real(dp), parameter :: M_pol = 2.0_dp, N_tor = nfp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.05_dp, eps_1 = -0.0375_dp
    type(anti_sigma_field_t) :: field

    integer :: this
    integer, parameter :: n_eta = 10
    real(dp), dimension(:), allocatable :: eta_grid

    logical, parameter :: should_plot = .false.
    logical :: test_failed

    test_failed = .false.

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)

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

    eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)

    do this = 1, n_fieldlines
        call calc_eta_integrand(field, &
                                fieldlines(this), &
                                eta_grid, &
                                eta_integrands(this))
        if (.not. size(eta_integrands(this)%eta_grid) .eq. n_eta) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_eta_integrand failed: eta grid not initialized"
            test_failed = .true.
        end if
        if (.not. size(eta_integrands(this)%F1) .eq. n_eta) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_eta_integrand failed: size(F1) =/= size(eta)"
            test_failed = .true.
        end if
        if (.not. size(eta_integrands(this)%F2) .eq. n_eta) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_eta_integrand failed: size(F2) =/= size(eta)"
            test_failed = .true.
        end if
        if (.not. size(eta_integrands(this)%F3) .eq. n_eta) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_eta_integrand failed: size(F3) =/= size(eta)"
            test_failed = .true.
        end if
    end do

    deallocate (eta_grid)

    if (should_plot) then
        call plot_eta_integrands(eta_integrands, fieldlines)
    end if

    if (test_failed) then
        error stop
    else
        if (.not. average_eta_integrands_sum_equal_F(field, &
                                                     fieldlines, &
                                                     M_pol, &
                                                     N_tor, &
                                                     eta_integrands)) then
            test_failed = .true.
        end if
    end if

    if (test_failed) error stop

end program test_calc_eta_integrand
