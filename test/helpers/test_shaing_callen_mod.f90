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
        use make_fieldline, only: make_flock_of_fieldlines
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

        call make_flock_of_fieldlines(flock, &
                                      xi_0, &
                                      iota, &
                                      circular_tok_field, &
                                      M_pol, &
                                      N_tor, &
                                      nfp)
        found = calc_trapped_fraction(circular_tok_field, flock, n_eta)

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
        trapped_fraction = calc_trapped_fraction(qs_field, qs_flock, n_eta)
        trapped_fraction_qs = calc_quasi_symmetric_trapped_fraction(qs_field, &
                                                                    eta_b, &
                                                                    n_eta)
        if (not_same(trapped_fraction, trapped_fraction_qs, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_trapped_fraction_against_qs failed:", &
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

    subroutine test_calc_avg_normalized_B_squared_dphimax_dxi0(qs_flock, &
                                                               test_failed)
        use shaing_callen_mod, only: calc_avg_normalized_B_squared_dphimax_dxi0
        type(flock_of_fieldlines_t), intent(in) :: qs_flock
        logical, intent(inout) :: test_failed

        real(dp), parameter :: reltol = 1e-12, abstol = 0.0_dp
        real(dp), dimension(size(qs_flock%fieldlines)) :: well_lengths
        real(dp) :: avg_B_squared, found, analytic

        real(dp) :: M_pol, nfp

        M_pol = qs_flock%M_pol
        nfp = qs_flock%nfp

        avg_B_squared = sum(qs_flock%fieldlines%phi_max(2) - &
                            qs_flock%fieldlines%phi_max(1))
        avg_B_squared = avg_B_squared/ &
                        sum(qs_flock%fieldlines%integral_one_over_B_squared)
        analytic = avg_B_squared*qs_flock%eta_b**2.0_dp*M_pol/nfp
        found = calc_avg_normalized_B_squared_dphimax_dxi0(qs_flock)

        if (not_same(found, analytic, reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_avg_normalized_B_squared_dphimax_dxi0 failed: ", &
                "for quasi-symmetric fields M_pol = ", M_pol, " nfp = ", nfp
            print *, "found = ", found
            print *, "analytic = ", analytic
            print *, "relative error = ", abs(1.0_dp - found/analytic)
            test_failed = .true.
        end if
    end subroutine test_calc_avg_normalized_B_squared_dphimax_dxi0

    subroutine test_calc_avg_normalized_lambda_dphimax_dxi0(qs_field, &
                                                            qs_flock, &
                                                            test_failed)
        use shaing_callen_mod, only: calc_avg_normalized_lambda_dphimax_dxi0
        use shaing_callen_mod, only: calc_avg_B_squared_over_avg_lambda
        use shaing_callen_integration, only: get_eta_integration_grid
        use shaing_callen_integration, only: integrate_over_eta_grid
        use fieldline_integrands, only: calc_lambda_squared
        use make_fieldline, only: get_global_B_max
        class(field_t), intent(in) :: qs_field
        type(flock_of_fieldlines_t), intent(in) :: qs_flock
        logical, intent(inout) :: test_failed

        integer, parameter :: n_eta = 10
        real(dp), parameter :: reltol = 1e-12, abstol = 0.0_dp
        real(dp), dimension(size(qs_flock%fieldlines)) :: well_lengths
        real(dp) :: avg_B_squared
        real(dp), dimension(n_eta) :: eta_grid, avg_lambda, lambda_max
        real(dp), dimension(n_eta) :: integrand
        real(dp), dimension(n_eta) :: found, analytic
        real(dp) :: found_integral, analytic_integral

        real(dp) :: M_pol, nfp, B_globalmax
        integer :: this

        M_pol = qs_flock%M_pol
        nfp = qs_flock%nfp
        B_globalmax = get_global_B_max(qs_flock%fieldlines)

        eta_grid = get_eta_integration_grid(qs_flock%eta_b, n_eta)
        avg_B_squared = sum(qs_flock%fieldlines%phi_max(2) - &
                            qs_flock%fieldlines%phi_max(1))
        avg_B_squared = avg_B_squared/ &
                        sum(qs_flock%fieldlines%integral_one_over_B_squared)
        avg_lambda = avg_B_squared/calc_avg_B_squared_over_avg_lambda(qs_field, &
                                                                      qs_flock, &
                                                                      eta_grid)
        do this = 1, n_eta
            lambda_max(this) = sqrt(calc_lambda_squared(B_globalmax, eta_grid(this)))
        end do
        analytic = avg_lambda/lambda_max*M_pol/nfp
        found = calc_avg_normalized_lambda_dphimax_dxi0(qs_field, &
                                                        qs_flock, &
                                                        eta_grid)

        if (not_same(found, analytic, reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_avg_normalized_lambda_dphimax_dxi0 failed: ", &
                "as function of eta"
            print *, "found = ", found
            print *, "analytic = ", analytic
            print *, "relative error = ", abs(1.0_dp - found/analytic)
            print *, "max relative error = ", maxval(abs(1.0_dp - found/analytic))
            test_failed = .true.
        end if

        integrand = eta_grid*analytic/avg_lambda
        analytic_integral = integrate_over_eta_grid(eta_grid, integrand)

        integrand = eta_grid*found/avg_lambda
        found_integral = integrate_over_eta_grid(eta_grid, integrand)

        if (not_same(found_integral, analytic_integral, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_avg_normalized_lambda_dphimax_dxi0 failed: ", &
                "integrals over eta"
            print *, "found = ", found_integral
            print *, "analytic = ", analytic_integral
            print *, "relative error = ", abs(1.0_dp - found_integral/analytic_integral)
            test_failed = .true.
        end if
    end subroutine test_calc_avg_normalized_lambda_dphimax_dxi0

    subroutine test_get_non_omnigenous_remainder(qs_field, qs_flock, test_failed)
        use shaing_callen_mod, only: get_non_omnigenous_remainder
        class(field_t), intent(in) :: qs_field
        type(flock_of_fieldlines_t), intent(in) :: qs_flock
        logical, intent(inout) :: test_failed

        real(dp), parameter :: reltol = 0.0_dp, const = 2.5e-1
        real(dp) :: abstol
        integer, parameter, dimension(5) :: n_etas = [50, 100, 200, 400, 800]

        real(dp) :: found, analytic
        integer :: this
        integer :: n_eta

        analytic = 0.0_dp
        do this = 1, size(n_etas)
            n_eta = n_etas(this)
            abstol = const/real(n_eta, kind=dp)**2.0_dp
            found = get_non_omnigenous_remainder(qs_field, qs_flock, n_eta)
            if (not_same(found, analytic, reltol_in=reltol, abstol_in=abstol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_get_non_omnigenous_remainder failed: quasi-symmetric"
                print *, "n_eta = ", n_eta
                print *, "found = ", found
                print *, "analytic = ", analytic
                print *, "abs error = ", abs(analytic - found)
                print *, "expected error = ", abstol
                test_failed = .true.
            end if
        end do
    end subroutine test_get_non_omnigenous_remainder

end module test_shaing_callen_mod
