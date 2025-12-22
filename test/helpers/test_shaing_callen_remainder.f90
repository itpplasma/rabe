module test_shaing_callen_remainder
    use constants, only: dp, pi
    use utils, only: not_same
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    use shaing_callen_integration, only: get_eta_integration_grid
    use test_shaing_callen_mod, only: calc_quasi_symmetric_trapped_fraction

    implicit none

contains

    subroutine test_avg_B_sq_antider_dBdtheta_over_B_cubed(qs_field, &
                                                           qs_fieldlines, &
                                                           test_failed)
     use shaing_callen_remainder, only: calc_avg_B_squared_antider_dBdtheta_over_B_cubed
        use shaing_callen_remainder, only: calc_avg_B_squared
        class(field_t), intent(in) :: qs_field
        type(fieldline_t), dimension(:), intent(in) :: qs_fieldlines
        logical, intent(inout) :: test_failed

        real(dp), parameter :: reltol = 1e-4, abstol = 0.0_dp

        real(dp) :: average, average_qs
        real(dp) :: avg_B_squared
        real(dp) :: M_pol, N_tor, iota

        average = calc_avg_B_squared_antider_dBdtheta_over_B_cubed(qs_field, &
                                                                   qs_fieldlines)
        avg_B_squared = calc_avg_B_squared(qs_fieldlines)
        average_qs = -0.5_dp*(1.0_dp - avg_B_squared*qs_fieldlines(1)%eta_b**2.0_dp)
        M_pol = qs_fieldlines(1)%M_pol
        N_tor = qs_fieldlines(1)%N_tor
        iota = qs_fieldlines(1)%iota
        average_qs = M_pol/(M_pol*iota - N_tor)*average_qs
        if (not_same(average, average_qs, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_avg_B_sq_antider_dBdtheta_over_B_cubed failed: ", &
                "<B^2 antider(dBdtheta/B^3)>"
            print *, "general = ", average
            print *, "quasi-symmetric = ", average_qs
            print *, "relative error = ", abs(1.0_dp - average/average_qs)
            test_failed = .true.
        end if
    end subroutine test_avg_B_sq_antider_dBdtheta_over_B_cubed

    subroutine test_trapped_fraction_prime_against_qs(qs_field, &
                                                      qs_fieldlines, &
                                                      test_failed)
        use shaing_callen_remainder, only: calc_trapped_fraction_prime
        class(field_t), intent(in) :: qs_field
        type(fieldline_t), dimension(:), intent(in) :: qs_fieldlines
        logical, intent(inout) :: test_failed

        real(dp), parameter :: reltol = 1e-3, abstol = 0.0_dp

        real(dp) :: eta_b
        integer, parameter :: n_eta = 100
        real(dp) :: trapped_fraction_prime, trapped_fraction_prime_qs
        real(dp) :: M_pol, N_tor, iota

        eta_b = qs_fieldlines(1)%eta_b
        trapped_fraction_prime = calc_trapped_fraction_prime(qs_field, &
                                                             qs_fieldlines, &
                                                             n_eta)
        trapped_fraction_prime_qs = calc_quasi_symmetric_trapped_fraction(qs_field, &
                                                                          eta_b, &
                                                                          n_eta)
        M_pol = qs_fieldlines(1)%M_pol
        N_tor = qs_fieldlines(1)%N_tor
        iota = qs_fieldlines(1)%iota
        trapped_fraction_prime_qs = M_pol/(M_pol*iota - N_tor)*trapped_fraction_prime_qs
        if (not_same(trapped_fraction_prime, trapped_fraction_prime_qs, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_trapped_fraction_prime_against_qs failed: ", &
                "trapped fraction prime expression"
            print *, "general = ", trapped_fraction_prime
            print *, "quasi-symmetric = ", trapped_fraction_prime_qs
            print *, "relative error = ", abs(1.0_dp - &
                                       trapped_fraction_prime/trapped_fraction_prime_qs)
            test_failed = .true.
        end if
    end subroutine test_trapped_fraction_prime_against_qs

    subroutine test_calc_avg_normalized_B_sq_dphimax_dxi0(qs_fieldlines, &
                                                          test_failed)
        use shaing_callen_remainder, only: calc_avg_normalized_B_squared_dphimax_dxi0
        type(fieldline_t), dimension(:), intent(in) :: qs_fieldlines
        logical, intent(inout) :: test_failed

        real(dp), parameter :: reltol = 1e-12, abstol = 0.0_dp
        real(dp) :: avg_B_squared, found, analytic

        real(dp) :: M_pol, N_tor

        M_pol = qs_fieldlines(1)%M_pol
        N_tor = qs_fieldlines(1)%N_tor

        avg_B_squared = sum(qs_fieldlines%phi_max(2) - qs_fieldlines%phi_max(1))
        avg_B_squared = avg_B_squared/sum(qs_fieldlines%integral_one_over_B_squared)
        analytic = avg_B_squared*qs_fieldlines(1)%eta_b**2.0_dp*M_pol/N_tor
        found = calc_avg_normalized_B_squared_dphimax_dxi0(qs_fieldlines)

        if (not_same(found, analytic, reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_calc_avg_normalized_B_sq_dphimax_dxi0 failed: ", &
                "for quasi-symmetric fields M_pol = ", M_pol, " N_tor = ", N_tor
            print *, "found = ", found
            print *, "analytic = ", analytic
            print *, "relative error = ", abs(1.0_dp - found/analytic)
            test_failed = .true.
        end if
    end subroutine test_calc_avg_normalized_B_sq_dphimax_dxi0

    subroutine test_calc_avg_normalized_lambda_dphimax_dxi0(qs_field, &
                                                            qs_fieldlines, &
                                                            test_failed)
        use shaing_callen_remainder, only: calc_avg_normalized_lambda_dphimax_dxi0
        use shaing_callen_mod, only: calc_avg_B_squared_over_avg_lambda
        use shaing_callen_integration, only: get_eta_integration_grid
        use shaing_callen_integration, only: integrate_over_eta_grid
        use fieldline_integrands, only: calc_lambda_squared
        class(field_t), intent(in) :: qs_field
        type(fieldline_t), dimension(:), intent(in) :: qs_fieldlines
        logical, intent(inout) :: test_failed

        integer, parameter :: n_eta = 10
        real(dp), parameter :: reltol = 1e-12, abstol = 0.0_dp
        real(dp), dimension(size(qs_fieldlines)) :: well_lengths
        real(dp) :: avg_B_squared
        real(dp), dimension(n_eta) :: eta_grid, avg_lambda, lambda_max
        real(dp), dimension(n_eta) :: integrand
        real(dp), dimension(n_eta) :: found, analytic
        real(dp) :: found_integral, analytic_integral

        real(dp) :: M_pol, nfp, B_globalmax
        integer :: this

        M_pol = qs_fieldlines(1)%M_pol
        nfp = qs_fieldlines(1)%nfp
        B_globalmax = 1.0_dp/qs_fieldlines(1)%eta_b

        eta_grid = get_eta_integration_grid(qs_fieldlines(1)%eta_b, n_eta)
        avg_B_squared = sum(qs_fieldlines%phi_max(2) - qs_fieldlines%phi_max(1))
        avg_B_squared = avg_B_squared/sum(qs_fieldlines%integral_one_over_B_squared)
        avg_lambda = avg_B_squared/calc_avg_B_squared_over_avg_lambda(qs_field, &
                                                                      qs_fieldlines, &
                                                                      eta_grid)
        do this = 1, n_eta
            lambda_max(this) = sqrt(calc_lambda_squared(B_globalmax, eta_grid(this)))
        end do
        analytic = avg_lambda/lambda_max*M_pol/nfp
        found = calc_avg_normalized_lambda_dphimax_dxi0(qs_field, &
                                                        qs_fieldlines, &
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

    subroutine test_get_non_omnigenous_remainder(qs_field, qs_fieldlines, test_failed)
        use shaing_callen_remainder, only: get_non_omnigenous_remainder
        class(field_t), intent(in) :: qs_field
        type(fieldline_t), dimension(:), intent(in) :: qs_fieldlines
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
            found = get_non_omnigenous_remainder(qs_field, qs_fieldlines, n_eta)
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

    subroutine test_get_non_omnigenous_remainders(qs_field, &
                                                  qs_fieldlines, &
                                                  test_failed)
        use shaing_callen_remainder, only: get_non_omnigenous_remainder_magnetic
        use shaing_callen_remainder, only: get_non_omnigenous_remainder_pitch
        class(field_t), intent(in) :: qs_field
        type(fieldline_t), dimension(:), intent(in) :: qs_fieldlines
        logical, intent(inout) :: test_failed

        real(dp), parameter :: reltol = 0.0_dp
        real(dp), parameter :: abstol_magnetics = 1e-16
        real(dp), parameter :: abstol_pitch = 1e-14

        integer, parameter :: n_eta = 100

        real(dp) :: found, analytic
        integer :: this

        analytic = 0.0_dp

        found = get_non_omnigenous_remainder_magnetic(qs_fieldlines)
        if (not_same(found, analytic, &
                     reltol_in=reltol, abstol_in=abstol_magnetics)) then
            print *, "-------------------------------------------------------------"
            print *, "test_get_non_omnigenous_remainders failed: magnetic part"
            print *, "found = ", found
            print *, "analytic = ", analytic
            print *, "abs error = ", abs(analytic - found)
            test_failed = .true.
        end if

        found = get_non_omnigenous_remainder_pitch(qs_field, &
                                                   qs_fieldlines, &
                                                   n_eta)
        if (not_same(found, analytic, &
                     reltol_in=reltol, abstol_in=abstol_pitch)) then
            print *, "-------------------------------------------------------------"
            print *, "test_get_non_omnigenous_remainders failed: pitch part"
            print *, "found = ", found
            print *, "analytic = ", analytic
            print *, "abs error = ", abs(analytic - found)
            test_failed = .true.
        end if
    end subroutine test_get_non_omnigenous_remainders

    subroutine test_limit_cancelation(fieldlines, test_failed)
        use shaing_callen_integration, only: get_eta_integration_grid
        use shaing_callen_integration, only: integrate_over_eta_grid
        use shaing_callen_remainder, only: calc_lambda_max
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        logical, intent(inout) :: test_failed

        real(dp), parameter :: reltol = 0.0_dp, const = 1.23_dp
        real(dp) :: abstol
        integer, parameter, dimension(5) :: n_etas = [50, 100, 200, 400, 800]

        integer :: n_eta
        real(dp) :: eta_b
        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: lambda_max
        real(dp) :: remainder, analytic

        integer :: this

        analytic = 0.0_dp
        eta_b = fieldlines(1)%eta_b
        do this = 1, size(n_etas)
            n_eta = n_etas(this)
            abstol = const/real(n_eta, kind=dp)**2.0_dp
            allocate (lambda_max(n_eta))
            eta_grid = get_eta_integration_grid(eta_b, n_eta)
            lambda_max = calc_lambda_max(eta_b, eta_grid)
            remainder = 0.75_dp*integrate_over_eta_grid(eta_grid, eta_grid/lambda_max)
            remainder = remainder - eta_b**2.0_dp
            if (not_same(remainder, analytic, reltol_in=reltol, abstol_in=abstol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_limit_cancelation failed:"
                print *, "n_eta = ", n_eta
                print *, "found = ", remainder
                print *, "analytic = ", analytic
                print *, "abs error = ", abs(analytic - remainder)
                print *, "expected error = ", abstol
                test_failed = .true.
            end if
            deallocate (lambda_max)
            deallocate (eta_grid)
        end do

    end subroutine test_limit_cancelation

end module test_shaing_callen_remainder
