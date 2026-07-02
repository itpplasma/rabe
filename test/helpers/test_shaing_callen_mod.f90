module test_shaing_callen_mod
    use constants, only: dp, pi
    use utils, only: not_same
    use fieldline_mod, only: flock_of_fieldlines_t
    use field_base, only: field_t
    use shaing_callen_integration, only: get_eta_integration_grid

    implicit none

contains

    subroutine test_trapped_fraction_against_circular_tokamak(test_failed)
        use shaing_callen_mod, only: calc_trapped_fraction
        use mock_field, only: mock_field_t
        use make_fieldline, only: make_flock_from_labels
        use utils, only: linspace

        logical, intent(inout) :: test_failed

        integer, parameter :: n_fieldlines = 50
        real(dp), dimension(n_fieldlines + 1) :: temp
        real(dp), dimension(n_fieldlines) :: xi_0
        type(mock_field_t) :: circular_tok_field
        type(flock_of_fieldlines_t) :: flock

        real(dp), parameter :: B_0 = 1.0_dp, eps = -0.001_dp
        real(dp), parameter :: M_pol = 1.0_dp, N_tor = 0.0_dp, nfp = 1.0_dp
        real(dp), parameter :: iota = 1.05_dp
        integer, parameter :: n_eta = 100

        real(dp), parameter :: reltol = 24.0_dp/real(n_eta, kind=dp)**2.0_dp
        real(dp), parameter :: abstol = 0.0_dp
        real(dp) :: found
        real(dp), parameter :: analytical_approx = 1.462_dp*sqrt(abs(eps))

        call circular_tok_field%mock_field_init(M_pol, N_tor, B_0, eps)

        call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
        xi_0 = temp(1:n_fieldlines)

        call make_flock_from_labels(flock, &
                                    xi_0, &
                                    iota, &
                                    circular_tok_field, &
                                    M_pol, &
                                    N_tor, &
                                    nfp)
        found = calc_trapped_fraction(flock, circular_tok_field, n_eta)

        if (not_same(found, analytical_approx, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_trapped_fraction_circular_tokamak failed:", &
                "trapped fraction expression"
            print *, "found = ", found
            print *, "approx analytic = ", analytical_approx
            print *, "relative error = ", abs(1.0_dp - found/analytical_approx)
            print *, "expected error = ", reltol
            test_failed = .true.
        end if
    end subroutine test_trapped_fraction_against_circular_tokamak

    subroutine test_trapped_fraction_against_qs(qs_field, qs_flock, test_failed)
        use shaing_callen_mod, only: calc_trapped_fraction
        class(field_t), intent(in) :: qs_field
        type(flock_of_fieldlines_t), intent(in) :: qs_flock
        logical, intent(inout) :: test_failed

        real(dp), parameter :: reltol = 1e-6, abstol = 0.0_dp

        real(dp) :: eta_b
        integer, parameter :: n_eta = 100
        real(dp) :: trapped_fraction, trapped_fraction_qs

        eta_b = qs_flock%eta_b
        trapped_fraction = calc_trapped_fraction(qs_flock, qs_field, n_eta)
        trapped_fraction_qs = calc_quasi_symmetric_trapped_fraction(qs_field, &
                                                                    eta_b, &
                                                                    n_eta)
        if (not_same(trapped_fraction, trapped_fraction_qs, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_trapped_fraction_against_qs failed: ", &
                "trapped fraction expression"
            print *, "general = ", trapped_fraction
            print *, "quasi-symmetric = ", trapped_fraction_qs
            print *, "relative error = ", abs(1.0_dp - &
                                              trapped_fraction/trapped_fraction_qs)
            test_failed = .true.
        end if
    end subroutine test_trapped_fraction_against_qs

    function calc_quasi_symmetric_trapped_fraction(field, &
                                                   eta_b, &
                                                   n_eta) result(trapped_fraction)
        use shaing_callen_integration, only: get_eta_integration_grid
        use shaing_callen_integration, only: integrate_over_eta_grid
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: eta_b
        integer, intent(in) :: n_eta
        real(dp) :: trapped_fraction

        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(n_eta) :: integrand
        real(dp) :: eta, average

        integer :: this

        eta_grid = get_eta_integration_grid(eta_b, n_eta)
        do this = 1, n_eta
            eta = eta_grid(this)
            average = get_theta_average_lambda_over_B_squared(field, eta)
            integrand(this) = eta/average
        end do
        trapped_fraction = 1.0_dp - 0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                    integrand)
    end function calc_quasi_symmetric_trapped_fraction

    function get_theta_average_lambda_over_B_squared(field, eta) result(average)
        use fieldline_integrands, only: lambda_over_B_squared
        use utils, only: linspace

        class(field_t), intent(in) :: field
        real(dp), intent(in) :: eta
        real(dp) :: average

        integer, parameter :: n_theta = 100
        real(dp), dimension(n_theta) :: thetas
        real(dp), dimension(n_theta - 1) :: integrand
        real(dp), parameter :: phi = 0.0_dp

        integer :: this

        call linspace(0.0_dp, 2.0_dp*pi, n_theta, thetas)
        do this = 1, n_theta - 1
            integrand(this) = lambda_over_B_squared(field, &
                                                    thetas(this), &
                                                    phi, &
                                                    eta)
        end do
        average = sum(integrand)/real(n_theta - 1, kind=dp)
    end function get_theta_average_lambda_over_B_squared

end module test_shaing_callen_mod
