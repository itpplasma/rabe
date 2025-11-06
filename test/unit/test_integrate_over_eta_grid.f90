program test_integrate_over_eta_grid
    use constants, only: dp
    use utils, only: not_same
    use shaing_callen_integration, only: get_eta_integration_grid
    use shaing_callen_integration, only: integrate_over_eta_grid

    implicit none

    real(dp), parameter :: eta_b = 2.0_dp
    real(dp), dimension(:), allocatable :: eta_grid, integrand
    real(dp), parameter :: integral = 4.0_dp*sqrt(eta_b)*(0.5_dp*log(eta_b) - 1.0_dp)
    real(dp) :: found_integral

    integer :: this, that

    real(dp), parameter :: abstol = 0.0_dp
    integer, parameter :: n_etas(6) = [100, 200, 400, 800, 1600, 3200]
    integer :: n_eta
    real(dp) :: reltol
    logical :: test_failed

    test_failed = .false.

    reltol = 6e-4
    do this = 1, size(n_etas)
        n_eta = n_etas(this)

        allocate (eta_grid(n_eta), integrand(n_eta))
        eta_grid = get_eta_integration_grid(eta_b, n_eta)
        do that = 1, n_eta
            integrand(that) = trial_func(eta_grid(that))
        end do
        found_integral = integrate_over_eta_grid(eta_grid, integrand)
        deallocate (eta_grid, integrand)
        if (not_same(integral, &
                     found_integral, &
                     reltol_in=reltol, &
                     abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_integrate_over_eta_grid failed: n_eta = ", n_eta
            print *, "found: ", found_integral
            print *, "analytic: ", integral
            print *, "relative error: ", abs(1.0_dp - found_integral/integral)
            print *, "expected error: ", reltol
            test_failed = .true.
        end if

        if (this < size(n_etas)) then
            reltol = reltol* &
                     real(n_etas(this), kind=dp)**2.0_dp/ &
                     real(n_etas(this + 1), kind=dp)**2.0_dp
        end if
    end do

    if (test_failed) error stop

contains

    function trial_func(eta) result(res)
        real(dp), intent(in) :: eta
        real(dp) :: res

        res = log(eta_b - eta)/sqrt(eta_b - eta)
    end function trial_func

end program test_integrate_over_eta_grid
