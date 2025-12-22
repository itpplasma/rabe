module shaing_callen_remainder
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    use shaing_callen_integration, only: get_eta_integration_grid
    use shaing_callen_integration, only: integrate_over_eta_grid
    implicit none

contains

    function calc_trapped_fraction_prime(field, &
                                         fieldlines, &
                                         n_eta) result(trapped_fraction_prime)
        use shaing_callen_mod, only: calc_avg_B_squared_over_avg_lambda
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        integer, intent(in) :: n_eta
        real(dp) :: trapped_fraction_prime

        integer :: this

        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda
        real(dp), dimension(:), allocatable :: F
        real(dp), dimension(:), allocatable :: integrand
        real(dp) :: avg_B_squared_antider_dBdtheta_over_B_cubed

        eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
        allocate (avg_B_squared_over_avg_lambda(n_eta))
        allocate (F(n_eta))
        allocate (integrand(n_eta))
        avg_B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(field, &
                                                                           fieldlines, &
                                                                           eta_grid)
        F = calc_F_prime(field, fieldlines, eta_grid)
        avg_B_squared_antider_dBdtheta_over_B_cubed = &
            calc_avg_B_squared_antider_dBdtheta_over_B_cubed(field, fieldlines)

        integrand = F*avg_B_squared_over_avg_lambda
        trapped_fraction_prime = -2.0_dp*avg_B_squared_antider_dBdtheta_over_B_cubed - &
                                 0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                 integrand)
        deallocate (integrand)
        deallocate (F)
        deallocate (avg_B_squared_over_avg_lambda)
        deallocate (eta_grid)

    end function calc_trapped_fraction_prime

    function calc_F_prime(field, fieldlines, eta_grid) result(F_prime)
        use shaing_callen_wrappers, only: this_field
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: F_prime

        integer :: this, n_fieldlines

        allocate (this_field, source=field)
        n_fieldlines = size(fieldlines)
        F_prime = 0.0_dp
        do this = 1, n_fieldlines
            F_prime = F_prime + calc_F_prime_for_fieldline(fieldlines(this), &
                                                           eta_grid)
        end do
        F_prime = F_prime/sum(fieldlines%integral_one_over_B_squared)
        deallocate (this_field)
    end function calc_F_prime

    function calc_F_prime_for_fieldline(fieldline, eta_grid) &
        result(F_prime)
        use shaing_callen_wrappers, only: wrapper_dBdtheta_over_lambda_cubed
        use shaing_callen_wrappers, only: wrapper_lambda_over_B_squared
        use shaing_callen_wrappers, only: this_eta, null_eta
        use shaing_callen_wrappers, only: this_fieldline, null_fieldline
        use fieldline_integrands, only: calc_lambda_squared
        use integrate, only: sum_trapez_1d
        use shaing_callen_integration, only: get_phi_integration_grid
        use shaing_callen_integration, only: integrate_over_phi_grid
        use shaing_callen_integration, only: cumintegrate_over_phi_grid

        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: F_prime

        real(dp), dimension(:), allocatable :: phi_grid
        real(dp) :: phi

        integer :: this, that, n_phi

        real(dp), dimension(:), allocatable :: antider_dBdtheta_over_lambda_cubed
        real(dp), dimension(:), allocatable :: phi_integrand_F

        this_fieldline = fieldline

        phi_grid = get_phi_integration_grid(fieldline)

        n_phi = size(phi_grid)
        allocate (antider_dBdtheta_over_lambda_cubed(n_phi))
        allocate (phi_integrand_F(n_phi))

        do this = 1, size(eta_grid)
            this_eta = eta_grid(this)
            antider_dBdtheta_over_lambda_cubed = cumintegrate_over_phi_grid(phi_grid, &
                                                     wrapper_dBdtheta_over_lambda_cubed)

            phi_integrand_F = 0.5_dp*this_eta**2.0_dp*antider_dBdtheta_over_lambda_cubed
            do that = 1, n_phi
                phi_integrand_F(that) = wrapper_lambda_over_B_squared(phi_grid(that)) &
                                        *phi_integrand_F(that)
            end do
            F_prime(this) = integrate_over_phi_grid(phi_grid, phi_integrand_F)
            this_eta = null_eta
        end do

        this_fieldline = null_fieldline
        deallocate (phi_grid)
        deallocate (antider_dBdtheta_over_lambda_cubed, phi_integrand_F)

    end function calc_F_prime_for_fieldline

    function calc_avg_B_squared_antider_dBdtheta_over_B_cubed(field, fieldlines) &
        result(res)
        use shaing_callen_wrappers, only: wrapper_dBdtheta_over_B_cubed
        use shaing_callen_wrappers, only: this_field
        use shaing_callen_wrappers, only: this_fieldline, null_fieldline
        use shaing_callen_integration, only: get_phi_integration_grid
        use shaing_callen_integration, only: cumintegrate_over_phi_grid
        use shaing_callen_integration, only: integrate_over_phi_grid
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: res

        integer :: this
        integer :: n_fieldlines
        real(dp), dimension(:), allocatable :: antider_dBdtheta_over_B_cubed
        real(dp), dimension(:), allocatable :: phi_grid

        n_fieldlines = size(fieldlines)
        allocate (this_field, source=field)
        res = 0.0_dp
        phi_grid = get_phi_integration_grid(fieldlines(1))
        allocate (antider_dBdtheta_over_B_cubed(size(phi_grid)))
        do this = 1, n_fieldlines
            this_fieldline = fieldlines(this)
            phi_grid = get_phi_integration_grid(this_fieldline)
            antider_dBdtheta_over_B_cubed = cumintegrate_over_phi_grid(phi_grid, &
                                                          wrapper_dBdtheta_over_B_cubed)
            res = res + integrate_over_phi_grid(phi_grid, antider_dBdtheta_over_B_cubed)
            this_fieldline = null_fieldline
        end do
        deallocate (antider_dBdtheta_over_B_cubed)
        deallocate (this_field)
        res = res/sum(fieldlines%integral_one_over_B_squared)
    end function calc_avg_B_squared_antider_dBdtheta_over_B_cubed

    function get_non_omnigenous_remainder(field, fieldlines, n_eta) result(remainder)
        use shaing_callen_integration, only: get_eta_integration_grid
        use shaing_callen_integration, only: integrate_over_eta_grid
        use shaing_callen_mod, only: calc_avg_B_squared_over_avg_lambda
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        integer, intent(in) :: n_eta
        real(dp) :: remainder

        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: integrand
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda

        eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
        allocate (integrand(n_eta))
        allocate (avg_B_squared_over_avg_lambda(n_eta))
        avg_B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(field, &
                                                                           fieldlines, &
                                                                           eta_grid)
        integrand = calc_avg_normalized_lambda_dphimax_dxi0(field, &
                                                            fieldlines, &
                                                            eta_grid)
        integrand = integrand*avg_B_squared_over_avg_lambda*eta_grid
        remainder = -(calc_avg_normalized_B_squared_dphimax_dxi0(fieldlines) - &
                      0.75_dp*integrate_over_eta_grid(eta_grid, integrand))
        deallocate (eta_grid)
        deallocate (integrand)
        deallocate (avg_B_squared_over_avg_lambda)

    end function get_non_omnigenous_remainder

    function get_non_omnigenous_remainder_magnetic(fieldlines) result(remainder)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: remainder

        real(dp) :: avg_B_squared
        real(dp) :: M_pol, nfp, B_max

        avg_B_squared = sum(fieldlines%phi_max(2) - fieldlines%phi_max(1))
        avg_B_squared = avg_B_squared/sum(fieldlines%integral_one_over_B_squared)

        M_pol = fieldlines(1)%M_pol
        nfp = fieldlines(1)%nfp
        B_max = 1.0_dp/fieldlines(1)%eta_b

        remainder = calc_avg_normalized_B_squared_dphimax_dxi0(fieldlines)
        remainder = remainder - avg_B_squared/B_max**2.0_dp*M_pol/nfp

    end function get_non_omnigenous_remainder_magnetic

    function get_non_omnigenous_remainder_pitch(field, fieldlines, n_eta) &
        result(remainder)
        use shaing_callen_integration, only: get_eta_integration_grid
        use shaing_callen_integration, only: integrate_over_eta_grid
        use shaing_callen_mod, only: calc_avg_B_squared_over_avg_lambda
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        integer, intent(in) :: n_eta
        real(dp) :: remainder

        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: integrand
        real(dp), dimension(:), allocatable :: lambda_max, omnigenous_integrand

        real(dp) :: avg_B_squared
        real(dp) :: M_pol, nfp, eta_b

        eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
        allocate (integrand(n_eta))
        allocate (lambda_max(n_eta), omnigenous_integrand(n_eta))

        integrand = calc_avg_B_squared_over_avg_lambda(field, fieldlines, eta_grid)
        integrand = integrand*calc_avg_normalized_lambda_dphimax_dxi0(field, &
                                                                      fieldlines, &
                                                                      eta_grid)
        integrand = 0.75_dp*integrand*eta_grid

        eta_b = fieldlines(1)%eta_b
        lambda_max = calc_lambda_max(eta_b, eta_grid)

        avg_B_squared = calc_avg_B_squared(fieldlines)

        M_pol = fieldlines(1)%M_pol
        nfp = fieldlines(1)%nfp
        omnigenous_integrand = 0.75_dp*avg_B_squared*eta_grid/lambda_max*M_pol/nfp
        remainder = integrate_over_eta_grid(eta_grid, integrand)
        remainder = remainder - integrate_over_eta_grid(eta_grid, omnigenous_integrand)

        deallocate (eta_grid)
        deallocate (integrand)
        deallocate (lambda_max)
        deallocate (omnigenous_integrand)

    end function get_non_omnigenous_remainder_pitch

    function calc_lambda_max(eta_b, eta_grid) result(lambda_max)
        real(dp), intent(in) :: eta_b
        real(dp), dimension(:) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: lambda_max

        lambda_max = 1.0_dp - eta_grid/eta_b
        if (any(lambda_max < 0.0_dp)) then
            print *, "Error in get_non_omnigenous_remainder_pitch:"
            print *, "1.0_dp - eta_grid/eta_b < 0"
            error stop
        end if
        lambda_max = sqrt(lambda_max)
    end function calc_lambda_max

    function calc_avg_B_squared(fieldlines) result(avg_B_squared)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: avg_B_squared

        avg_B_squared = sum(fieldlines%phi_max(2) - fieldlines%phi_max(1))
        avg_B_squared = avg_B_squared/sum(fieldlines%integral_one_over_B_squared)

    end function calc_avg_B_squared

    function calc_avg_normalized_B_squared_dphimax_dxi0(fieldlines) &
        result(res)
        use fourier, only: real_ft
        use fieldline_integrals, only: allocate_modes
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: res

        real(dp), dimension(size(fieldlines)) :: w_l, phi_l, xi_0
        real(dp), dimension(size(fieldlines)) :: dphimax_dxi0, well_lengths
        real(dp) :: M_pol, nfp

        M_pol = fieldlines(1)%M_pol
        nfp = fieldlines(1)%nfp
        phi_l = fieldlines%phi_max(1)
        xi_0 = fieldlines%xi_0
        w_l = phi_l - M_pol/nfp*xi_0
        dphimax_dxi0 = M_pol/nfp + calc_periodic_dydx(xi_0, w_l)
        well_lengths = fieldlines%phi_max(2) - fieldlines%phi_max(1)
        res = sum(well_lengths*dphimax_dxi0/fieldlines%B_max(1)**2.0_dp)
        res = res/sum(fieldlines%integral_one_over_B_squared)
    end function calc_avg_normalized_B_squared_dphimax_dxi0

    function calc_avg_normalized_lambda_dphimax_dxi0(field, fieldlines, eta_grid) &
        result(res)
        use fourier, only: real_ft
        use fieldline_integrals, only: allocate_modes
        use shaing_callen_wrappers, only: this_field
        use fieldline_integrands, only: calc_lambda_squared
        use shaing_callen_mod, only: calc_avg_lambda_over_B_squared
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: res

        real(dp), dimension(size(fieldlines)) :: w_l, phi_l, xi_0, B_l
        real(dp), dimension(size(fieldlines)) :: dphimax_dxi0
        real(dp), dimension(size(eta_grid)) :: lambda_l
        real(dp) :: M_pol, nfp
        integer :: this, that

        M_pol = fieldlines(1)%M_pol
        nfp = fieldlines(1)%nfp
        phi_l = fieldlines%phi_max(1)
        xi_0 = fieldlines%xi_0
        w_l = phi_l - M_pol/nfp*xi_0
        dphimax_dxi0 = M_pol/nfp + calc_periodic_dydx(xi_0, w_l)
        B_l = fieldlines%B_max(1)
        allocate (this_field, source=field)
        res = 0.0_dp
        do this = 1, size(fieldlines)
            do that = 1, size(eta_grid)
                lambda_l(that) = sqrt(calc_lambda_squared(B_l(this), eta_grid(that)))
            end do
            res = res + &
                  calc_avg_lambda_over_B_squared(fieldlines(this), eta_grid)* &
                  dphimax_dxi0(this)/lambda_l
        end do
        deallocate (this_field)
        res = res/sum(fieldlines%integral_one_over_B_squared)
    end function calc_avg_normalized_lambda_dphimax_dxi0

    function calc_periodic_dydx(x, y) result(dydx)
        use fourier, only: check_is_equidistant, check_has_correct_endpoints
        real(dp), dimension(:), intent(in) :: x, y
        real(dp), dimension(:), allocatable :: dydx

        real(dp), dimension(:), allocatable :: y_periodic
        real(dp) :: dx

        integer :: n

        call check_is_equidistant(x)
        call check_has_correct_endpoints(x)
        n = size(x)
        if (n < 3) then
            print *, "calc_periodic_dydx needs n > 3 points!"
            print *, "provided n: ", n
            error stop
        end if

        dx = x(2) - x(1)

        allocate (y_periodic(n + 2))
        y_periodic(1) = y(n)
        y_periodic(2:n + 1) = y
        y_periodic(n + 2) = y(1)

        allocate (dydx(n))
        dydx = 0.5_dp*(y_periodic(3:n + 2) - y_periodic(1:n))/dx
        deallocate (y_periodic)

    end function calc_periodic_dydx

end module shaing_callen_remainder
