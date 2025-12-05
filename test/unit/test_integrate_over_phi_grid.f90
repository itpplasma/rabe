program test_integrate_over_phi_grid
    use constants, only: dp, pi
    use utils, only: not_same
    use shaing_callen_integration, only: get_phi_integration_grid
    use shaing_callen_integration, only: integrate_over_phi_grid
    use fieldline_mod, only: fieldline_t

    implicit none

    type(fieldline_t) :: fieldline
    real(dp), dimension(2), parameter :: phi_limits = [0.0_dp, pi]
    real(dp), dimension(:), allocatable :: phi_grid, integrand
    real(dp), parameter :: integral = sqrt(8.0_dp)
    real(dp) :: found_integral

    integer :: this, that

    real(dp), parameter :: abstol = 0.0_dp
    integer, parameter :: n_phis(6) = [100, 200, 400, 800, 1600, 3200]
    integer :: n_phi
    real(dp) :: reltol
    logical :: test_failed

    test_failed = .false.
    fieldline%phi_max = phi_limits

    reltol = 6.6e-5
    do this = 1, size(n_phis)
        n_phi = n_phis(this)

        allocate (phi_grid(n_phi), integrand(n_phi))
        phi_grid = get_phi_integration_grid(fieldline, n_phi)
        do that = 1, n_phi
            integrand(that) = trial_func(phi_grid(that))
        end do
        found_integral = integrate_over_phi_grid(phi_grid, integrand)
        deallocate (phi_grid, integrand)
        if (not_same(integral, &
                     found_integral, &
                     reltol_in=reltol, &
                     abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_integrate_over_phi_grid failed: n_phi = ", n_phi
            print *, "found: ", found_integral
            print *, "analytic: ", integral
            print *, "relative error: ", abs(1.0_dp - found_integral/integral)
            print *, "expected error: ", reltol
            test_failed = .true.
        end if

        if (this < size(n_phis)) then
            reltol = reltol* &
                     real(n_phis(this), kind=dp)**2.0_dp/ &
                     real(n_phis(this + 1), kind=dp)**2.0_dp
        end if
    end do

    if (test_failed) error stop

contains

    function trial_func(phi) result(res)
        real(dp), intent(in) :: phi
        real(dp) :: res

        res = sqrt(1.0_dp - cos(phi))
    end function trial_func

end program test_integrate_over_phi_grid
