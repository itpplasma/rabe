module shaing_callen_mod
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    use integrate, only: integrate_1d_substituted, sum_trapez_1d
    use shaing_callen_integration, only: get_phi_integration_grid
    use shaing_callen_integration, only: integrate_over_phi_grid
    use shaing_callen_integration, only: cumint
    use shaing_callen_integration, only: integrate_over_eta_grid
    use shaing_callen_wrappers, only: wrapper_one_over_lambda
    use shaing_callen_wrappers, only: wrapper_lambda_over_B_squared
    implicit none

    type :: eta_integrand_t
        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: F1, F2, F3
    end type eta_integrand_t

    type :: shaing_callen_t
        real(dp) :: trapped_fraction
        real(dp) :: modified_trapped_fraction
    end type shaing_callen_t

contains

    function calc_shaing_callen(field, &
                                fieldlines, &
                                eta_grid) result(shaing_callen)
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:), intent(in) :: eta_grid
        type(shaing_callen_t) :: shaing_callen

        integer :: this

        integer :: n_eta
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda
        real(dp), dimension(:), allocatable :: F
        real(dp), dimension(:), allocatable :: integrand
        real(dp) :: trapped_fraction, modified_trapped_fraction
        real(dp) :: M_pol, N_tor, iota

        n_eta = size(eta_grid)
        allocate (avg_B_squared_over_avg_lambda(n_eta))
        allocate (F(n_eta))
        allocate (integrand(n_eta))
        avg_B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(field, &
                                                                           fieldlines, &
                                                                           eta_grid)
        F = calc_F(field, fieldlines, eta_grid)

        integrand = eta_grid*avg_B_squared_over_avg_lambda
        trapped_fraction = 1.0_dp - 0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                    integrand)

        integrand = F*avg_B_squared_over_avg_lambda
        M_pol = fieldlines(1)%M_pol
        N_tor = fieldlines(1)%N_tor
        iota = fieldlines(1)%iota
        modified_trapped_fraction = M_pol/(M_pol*iota - N_tor) - &
                                    0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                    integrand)
        deallocate (integrand)
        deallocate (F)
        deallocate (avg_B_squared_over_avg_lambda)
        shaing_callen%trapped_fraction = trapped_fraction
        shaing_callen%modified_trapped_fraction = modified_trapped_fraction

    end function calc_shaing_callen

    function calc_F(field, fieldlines, eta_grid) result(F)
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:) :: eta_grid
        real(dp), dimension(:), allocatable :: F

        integer :: this
        integer :: n_eta, n_fieldlines
        type(eta_integrand_t) :: eta_integrand
        real(dp) :: prefactor

        prefactor = fieldlines(1)%M_pol/ &
                    (fieldlines(1)%M_pol*fieldlines(1)%iota - fieldlines(1)%N_tor)

        n_eta = size(eta_grid)
        n_fieldlines = size(fieldlines)
        allocate (F(n_eta))
        F = 0.0_dp
        do this = 1, n_fieldlines
            call calc_eta_integrand(field, fieldlines(this), eta_grid, eta_integrand)
            F = F + &
                eta_integrand%F1*prefactor + &
                eta_integrand%F2 + &
                eta_integrand%F3
        end do
        F = F/sum(fieldlines%integral_one_over_B_squared)
    end function calc_F

    subroutine calc_eta_integrand(field, fieldline, eta_grid, eta_integrand)
        use fieldline_integrands, only: calc_lambda_squared
        use shaing_callen_wrappers, only: this_fieldline, this_eta, this_field
        use shaing_callen_wrappers, only: wrapper_lambda_dBdtheta_over_B_cubed
        use shaing_callen_wrappers, only: wrapper_dBdtheta_over_B_squared_lambda
        use shaing_callen_integration, only: cumintegrate_over_phi_grid

        type(fieldline_t), intent(in) :: fieldline
        class(field_t), intent(in) :: field
        real(dp), dimension(:), intent(in) :: eta_grid
        type(eta_integrand_t), intent(out) :: eta_integrand

        real(dp), dimension(:), allocatable :: phi_grid

        integer :: this, that, n_phi

        real(dp), dimension(:), allocatable :: antider_one_over_lambda
        real(dp), dimension(:), allocatable :: phi_integrand_F2
        real(dp), dimension(:), allocatable :: phi_integrand_F3
        real(dp) :: phi

        allocate (this_field, source=field)
        this_fieldline = fieldline

        phi_grid = get_phi_integration_grid(fieldline)

        call allocate_eta_integrands(eta_grid, eta_integrand)
        n_phi = size(phi_grid)
        allocate (antider_one_over_lambda(n_phi))
        allocate (phi_integrand_F2(n_phi))
        allocate (phi_integrand_F3(n_phi))

        do this = 1, size(eta_grid)
            this_eta = eta_grid(this)
            antider_one_over_lambda = cumintegrate_over_phi_grid(phi_grid, &
                                                                wrapper_one_over_lambda)
            do that = 1, n_phi
                phi = phi_grid(that)
                phi_integrand_F2(that) = wrapper_lambda_dBdtheta_over_B_cubed(phi)
                phi_integrand_F3(that) = wrapper_dBdtheta_over_B_squared_lambda(phi)
            end do
            phi_integrand_F2 = phi_integrand_F2*antider_one_over_lambda
            phi_integrand_F3 = phi_integrand_F3*antider_one_over_lambda
            eta_integrand%F1(this) = wrapper_lambda_over_B_squared(fieldline%phi_max(2))
          eta_integrand%F1(this) = eta_integrand%F1(this)*antider_one_over_lambda(n_phi)
            eta_integrand%F2(this) = integrate_over_phi_grid(phi_grid, phi_integrand_F2)
            eta_integrand%F3(this) = integrate_over_phi_grid(phi_grid, phi_integrand_F3)
        end do
        eta_integrand%F1 = eta_grid*eta_integrand%F1
        eta_integrand%F2 = 2.0_dp*eta_grid*eta_integrand%F2
        eta_integrand%F3 = 0.5_dp*eta_grid**2.0_dp*eta_integrand%F3

        deallocate (this_field)
        deallocate (phi_grid)

    end subroutine calc_eta_integrand

    subroutine allocate_eta_integrands(eta_grid, eta_integrand)
        real(dp), dimension(:), intent(in) :: eta_grid
        type(eta_integrand_t), intent(inout) :: eta_integrand

        integer :: n_eta

        if (allocated(eta_integrand%eta_grid)) deallocate (eta_integrand%eta_grid)
        if (allocated(eta_integrand%F1)) deallocate (eta_integrand%F1)
        if (allocated(eta_integrand%F2)) deallocate (eta_integrand%F2)
        if (allocated(eta_integrand%F3)) deallocate (eta_integrand%F3)

        n_eta = size(eta_grid)
        allocate (eta_integrand%eta_grid, source=eta_grid)
        allocate (eta_integrand%F1(n_eta))
        allocate (eta_integrand%F2(n_eta))
        allocate (eta_integrand%F3(n_eta))

    end subroutine allocate_eta_integrands

    function calc_alternative_F(field, fieldlines, eta_grid) result(F)
        use shaing_callen_wrappers, only: this_field
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: F

        integer :: this, n_fieldlines

        allocate (this_field, source=field)
        n_fieldlines = size(fieldlines)
        F = 0.0_dp
        do this = 1, n_fieldlines
            F = F + calc_alternative_F_for_fieldline(fieldlines(this), &
                                                     eta_grid)
        end do
        F = F/sum(fieldlines%integral_one_over_B_squared)
        deallocate (this_field)
    end function calc_alternative_F

    function calc_alternative_F_for_fieldline(fieldline, eta_grid) &
        result(F)
        use shaing_callen_wrappers, only: this_eta, null_eta
        use shaing_callen_wrappers, only: this_fieldline, null_fieldline
        use shaing_callen_integration, only: integral_dBdtheta_over_lambda_cubed
        use fieldline_integrands, only: calc_lambda_squared
        use integrate, only: sum_trapez_1d

        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: F

        real(dp) :: M_pol, N_tor
        real(dp), dimension(:), allocatable :: phi_grid
        real(dp) :: phi

        integer :: this, that, n_phi

        real(dp), dimension(:), allocatable :: antider_dBdtheta_over_lambda_cubed
        real(dp), dimension(:), allocatable :: phi_integrand_F

        this_fieldline = fieldline
        M_pol = fieldline%M_pol
        N_tor = fieldline%N_tor

        phi_grid = get_phi_integration_grid(fieldline)

        n_phi = size(phi_grid)
        allocate (antider_dBdtheta_over_lambda_cubed(n_phi))
        allocate (phi_integrand_F(n_phi))

        do this = 1, size(eta_grid)
            this_eta = eta_grid(this)
            antider_dBdtheta_over_lambda_cubed = cumint(phi_grid, &
                                                    integral_dBdtheta_over_lambda_cubed)

            phi_integrand_F = wrapper_one_over_lambda(fieldline%phi_max(1)) &
                              *M_pol/(N_tor - M_pol*fieldline%iota) &
                              - 0.5_dp*this_eta*antider_dBdtheta_over_lambda_cubed
            do that = 1, n_phi
                phi_integrand_F(that) = wrapper_lambda_over_B_squared(phi_grid(that)) &
                                        *phi_integrand_F(that)
            end do
            F(this) = sum_trapez_1d(phi_grid, phi_integrand_F)
            this_eta = null_eta
        end do
        F = -eta_grid*F

        this_fieldline = null_fieldline
        deallocate (phi_grid)
        deallocate (antider_dBdtheta_over_lambda_cubed, phi_integrand_F)

    end function calc_alternative_F_for_fieldline

    function calc_avg_B_squared_over_avg_lambda(field, fieldlines, eta_grid) &
        result(avg_B_squared_over_avg_lambda)
        use shaing_callen_wrappers, only: this_field
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda

        integer :: this
        integer :: n_eta, n_fieldlines
        real(dp) :: avg_well_length
        real(dp), dimension(:), allocatable :: avg_lambda_over_B_squared

        n_eta = size(eta_grid)
        n_fieldlines = size(fieldlines)
        allocate (avg_lambda_over_B_squared(n_eta))
        avg_well_length = 0.0_dp
        avg_lambda_over_B_squared = 0.0_dp
        allocate (this_field, source=field)
        do this = 1, n_fieldlines
            avg_lambda_over_B_squared = avg_lambda_over_B_squared + &
                              calc_avg_lambda_over_B_squared(fieldlines(this), eta_grid)
            avg_well_length = avg_well_length + &
                              fieldlines(this)%phi_max(2) - &
                              fieldlines(this)%phi_max(1)
        end do
        deallocate (this_field)
        allocate (avg_B_squared_over_avg_lambda(n_eta))

        avg_B_squared_over_avg_lambda = avg_well_length/ &
                                        avg_lambda_over_B_squared

        deallocate (avg_lambda_over_B_squared)

    end function calc_avg_B_squared_over_avg_lambda

    function calc_avg_lambda_over_B_squared(fieldline, eta_grid) &
        result(avg_lambda_over_B_squared)
        use shaing_callen_wrappers, only: this_eta, null_eta
        use shaing_callen_wrappers, only: this_fieldline, null_fieldline
        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: avg_lambda_over_B_squared

        integer :: this

        this_fieldline = fieldline

        do this = 1, size(eta_grid)
            this_eta = eta_grid(this)
            call integrate_1d_substituted(wrapper_lambda_over_B_squared, &
                                          fieldline%phi_max(1), &
                                          fieldline%phi_max(2), &
                                          avg_lambda_over_B_squared(this))
            this_eta = null_eta
        end do

        this_fieldline = null_fieldline

    end function calc_avg_lambda_over_B_squared

    function calc_shaing_callen_prime(field, &
                                      fieldlines, &
                                      eta_grid) result(shaing_callen)
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:), intent(in) :: eta_grid
        type(shaing_callen_t) :: shaing_callen

        integer :: this

        integer :: n_eta
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda
        real(dp), dimension(:), allocatable :: F_prime
        real(dp), dimension(:), allocatable :: integrand
        real(dp) :: trapped_fraction, modified_trapped_fraction
        real(dp) :: temp

        n_eta = size(eta_grid)
        allocate (avg_B_squared_over_avg_lambda(n_eta))
        allocate (integrand(n_eta))
        avg_B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(field, &
                                                                           fieldlines, &
                                                                           eta_grid)

        integrand = eta_grid*avg_B_squared_over_avg_lambda
        trapped_fraction = 1.0_dp - 0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                    integrand)

        allocate (F_prime(n_eta))

        F_prime = calc_F_prime(field, fieldlines, eta_grid)
        integrand = F_prime*avg_B_squared_over_avg_lambda
        modified_trapped_fraction = -2.0_dp*calc_avg_B_squared_antider_dBdtheta_over_B_cubed(field, fieldlines) &
                                    - 0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                      integrand)
      temp = -2.0_dp*calc_avg_B_squared_antider_dBdtheta_over_B_cubed(field, fieldlines)
        temp = temp*(fieldlines(1)%iota - fieldlines(1)%N_tor/fieldlines(1)%M_pol)
        temp = temp + fieldlines(1)%eta_b**2.0_dp* &
               (fieldlines(1)%phi_max(2) - fieldlines(1)%phi_max(1))/ &
               fieldlines(1)%integral_one_over_B_squared
        print *, temp
        deallocate (integrand)
        deallocate (F_prime)
        deallocate (avg_B_squared_over_avg_lambda)
        shaing_callen%trapped_fraction = trapped_fraction
        shaing_callen%modified_trapped_fraction = modified_trapped_fraction

    end function calc_shaing_callen_prime

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
        use shaing_callen_integration, only: integral_dBdtheta_over_lambda_cubed
        use shaing_callen_wrappers, only: this_eta, null_eta
        use shaing_callen_wrappers, only: this_fieldline, null_fieldline
        use fieldline_integrands, only: calc_lambda_squared
        use integrate, only: sum_trapez_1d

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
            antider_dBdtheta_over_lambda_cubed = cumint(phi_grid, &
                                                    integral_dBdtheta_over_lambda_cubed)

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
        use shaing_callen_wrappers, only: this_field
        use shaing_callen_wrappers, only: this_fieldline, null_fieldline
        use shaing_callen_integration, only: integral_dBdtheta_over_B_cubed
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
            antider_dBdtheta_over_B_cubed = cumint(phi_grid, &
                                                   integral_dBdtheta_over_B_cubed)
            res = res + integrate_over_phi_grid(phi_grid, antider_dBdtheta_over_B_cubed)
            this_fieldline = null_fieldline
        end do
        deallocate (antider_dBdtheta_over_B_cubed)
        deallocate (this_field)
        res = res/sum(fieldlines%integral_one_over_B_squared)

    end function calc_avg_B_squared_antider_dBdtheta_over_B_cubed

    function get_non_omnigenous_remainder(field, fieldlines, eta_grid) result(remainder)
        use shaing_callen_integration, only: integrate_over_eta_grid
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:) :: eta_grid
        real(dp) :: remainder

        real(dp), dimension(size(eta_grid)) :: integrand, B_squared_over_avg_lambda

        B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(field, &
                                                                       fieldlines, &
                                                                       eta_grid)
        integrand = calc_avg_normalized_lambda_dphimax_dxi0(field, &
                                                            fieldlines, &
                                                            eta_grid)
        integrand = integrand*B_squared_over_avg_lambda*eta_grid
        remainder = calc_avg_normalized_B_squared_dphimax_dxi0(fieldlines) - &
                    0.75_dp*integrate_over_eta_grid(eta_grid, integrand)

    end function get_non_omnigenous_remainder

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
        real(dp), dimension(:), intent(in) :: x, y
        real(dp), dimension(:), allocatable :: dydx

        real(dp), dimension(:), allocatable :: y_periodic
        real(dp) :: dx

        integer :: n

        call check_is_equidistant(x)
        call check_has_periodic_endpoints(x)
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

    subroutine check_is_equidistant(x)
        real(dp), dimension(:), intent(in) :: x
        logical :: is_equidistant

        real(dp), parameter :: retol = 1e-13, abstol = 1e-14
        real(dp), dimension(size(x) - 1) :: dx
        real(dp) :: tol
        integer :: N

        N = size(x)

        dx = x(2:N) - x(1:N - 1)
        is_equidistant = all(abs(dx - dx(1)) < retol*dx(1) + abstol)

        if (.not. is_equidistant) then
            print *, "Input x has to be equidistant!"
            print *, "violation by ", maxval(abs(dx - dx(1)))
            error stop
        end if
    end subroutine check_is_equidistant

    subroutine check_has_periodic_endpoints(x)
        use constants, only: pi
        real(dp), dimension(:), intent(in) :: x
        logical :: has_correct_range

        real(dp), parameter :: tol = 1e-13
        real(dp) :: correct_range, range
        integer :: N

        N = size(x)
        range = x(N) - x(1)
        correct_range = 2.0_dp*pi*(1 - 1/real(N, kind=dp))
        has_correct_range = abs(correct_range - range) < tol*correct_range
        if (.not. has_correct_range) then
            print *, "Input x has wrong endpoints!"
            print *, "actual: ", x(1), x(N)
            print *, "required: ", x(1), x(1) + correct_range
            error stop
        end if
    end subroutine check_has_periodic_endpoints

end module shaing_callen_mod
