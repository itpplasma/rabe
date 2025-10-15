module shaing_callen_mod
    use constants, only: dp, pi
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    use integrate, only: integrate_1d_substituted, sum_trapez_1d
    implicit none

    type :: eta_integrand_t
        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: F1, F2, F3
        real(dp), dimension(:), allocatable :: integral_lambda_over_B_squared
    end type eta_integrand_t

    abstract interface
        function defined_integral_i(x_start, x_end)
            use constants, only: dp
            real(dp), intent(in) :: x_start, x_end
            real(dp) :: defined_integral_i
        end function defined_integral_i
    end interface

    class(field_t), allocatable, private :: this_field
    type(fieldline_t), private :: this_fieldline
    real(dp), private :: this_eta

contains

    function calc_shaing_callen(fieldlines, eta_integrands) result(shaing_callen)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(eta_integrand_t), dimension(:), intent(in) :: eta_integrands
        real(dp) :: shaing_callen

        integer :: this

        integer :: n_eta, n_fieldlines
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda
        real(dp) :: trapped_fraction

        n_eta = size(eta_integrands(1)%eta_grid)
        n_fieldlines = size(fieldlines)
        allocate (avg_B_squared_over_avg_lambda(n_eta))
        avg_B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(fieldlines, &
                                                                         eta_integrands)
        deallocate (avg_B_squared_over_avg_lambda)

    end function calc_shaing_callen

    function calc_trapped_particle_fraction(fieldlines, &
                                            eta_integrands) result(trapped_fraction)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(eta_integrand_t), dimension(:), intent(in) :: eta_integrands
        real(dp) :: trapped_fraction

        integer :: this
        integer :: n_eta, n_fieldlines
        real(dp) :: avg_well_length
        real(dp), dimension(:), allocatable :: avg_lambda_over_B_squared
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda
        real(dp) :: dxi_0

    end function calc_trapped_particle_fraction

    function calc_avg_B_squared_over_avg_lambda(fieldlines, &
                                                eta_integrands) &
        result(avg_B_squared_over_av_lambda)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(eta_integrand_t), dimension(:), intent(in) :: eta_integrands
        real(dp) :: avg_B_squared_over_av_lambda

        integer :: this
        integer :: n_eta, n_fieldlines
        real(dp) :: avg_well_length
        real(dp), dimension(:), allocatable :: avg_lambda_over_B_squared
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda
        real(dp) :: dxi_0

        n_eta = size(eta_integrands(1)%eta_grid)
        n_fieldlines = size(fieldlines)
        allocate (avg_lambda_over_B_squared(n_eta))
        avg_well_length = 0.0_dp
        avg_lambda_over_B_squared = 0.0_dp
        do this = 1, n_fieldlines
            avg_lambda_over_B_squared = avg_lambda_over_B_squared + &
                                     eta_integrands(this)%integral_lambda_over_B_squared
            avg_well_length = avg_well_length + &
                              fieldlines(this)%phi_max(2) - &
                              fieldlines(this)%phi_max(1)
        end do
        dxi_0 = (fieldlines(n_fieldlines)%xi_0 - fieldlines(1)%xi_0)/ &
                real(n_fieldlines, kind=dp)
        avg_lambda_over_B_squared = avg_lambda_over_B_squared*dxi_0
        avg_well_length = avg_well_length*dxi_0
        allocate (avg_B_squared_over_avg_lambda(n_eta))

        avg_B_squared_over_avg_lambda = avg_well_length/ &
                                        avg_lambda_over_B_squared

        deallocate (avg_lambda_over_B_squared)
        deallocate (avg_B_squared_over_avg_lambda)

    end function calc_avg_B_squared_over_avg_lambda

    subroutine calc_eta_integrand(field, fieldline, eta_grid, eta_integrand)
        use fieldline_integrands, only: calc_lambda_squared

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
            call integrate_1d_substituted(wrapper_lambda_over_B_squared, &
                                          this_fieldline%phi_max(1), &
                                          this_fieldline%phi_max(2), &
                                     eta_integrand%integral_lambda_over_B_squared(this))
            antider_one_over_lambda = cumint(phi_grid, &
                                             integral_one_over_lambda)
            do that = 1, n_phi
                phi = phi_grid(that)
                phi_integrand_F2(that) = wrapper_lambda_dBdtheta_over_B_cubed(phi)
                phi_integrand_F3(that) = wrapper_dBdtheta_over_B_squared_lambda(phi)
            end do
            phi_integrand_F2 = phi_integrand_F2*antider_one_over_lambda
            phi_integrand_F3 = phi_integrand_F3*antider_one_over_lambda
            eta_integrand%F1(this) = wrapper_lambda_over_B_squared(fieldline%phi_max(2))
          eta_integrand%F1(this) = eta_integrand%F1(this)*antider_one_over_lambda(n_phi)
            eta_integrand%F2(this) = sum_trapez_1d(phi_grid, phi_integrand_F2)
            eta_integrand%F3(this) = sum_trapez_1d(phi_grid, phi_integrand_F3)
        end do
        eta_integrand%F1 = eta_grid*eta_integrand%F1
        eta_integrand%F2 = 2.0_dp*eta_grid*eta_integrand%F2
        eta_integrand%F3 = 0.5_dp*eta_grid**2.0_dp*eta_integrand%F3

        deallocate (this_field)
        deallocate (phi_grid)

    end subroutine calc_eta_integrand

    function get_eta_integration_grid(eta_b, n_eta) result(eta_grid)
        use utils, only: linspace
        real(dp), intent(in) :: eta_b
        integer, intent(in) :: n_eta
        real(dp), dimension(:), allocatable :: eta_grid

        real(dp), dimension(:), allocatable :: t

        allocate (eta_grid(n_eta), t(n_eta + 1))
        ! transform integral to variable t -> eta = eta_b - t**4
        ! In this transformed integral, the here considered integrands become
        ! zero at the upper limit eta = eta_b (t=0), due to the differential
        ! dt/deta = 4t**3 killing the logarithmic divergences. We therefore
        ! exlude eta = eta_b, t = 0, from the integral grid to avoid numerical
        ! instabilities due to the divergence.
        call linspace(eta_b**0.25_dp, 0.0_dp, n_eta + 1, t)
        eta_grid = eta_b - t(1:n_eta)**4.0_dp
        deallocate (t)
    end function get_eta_integration_grid

    function integrate_over_eta_grid(eta_grid, integrand) result(integral)
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(:), intent(in) :: integrand
        real(dp) :: integral

        integer :: n_eta
        real(dp), dimension(:), allocatable :: d_eta_dt
        real(dp) :: eta_b, dt

        n_eta = size(eta_grid)
        ! As last eta value = eta_b (1 - 1/n_eta**4) according to above grid
        eta_b = n_eta**4.0_dp/(n_eta**4.0_dp - 1.0_dp)*eta_grid(n_eta)
        allocate (d_eta_dt(n_eta))
        d_eta_dt = 4.0_dp*(eta_b - eta_grid)**0.75_dp
        dt = eta_b**0.25_dp/n_eta
        ! The considered integrands vanish at the upper limit (t=0) due to the
        ! differential 4t**3 (see comments in get_eta_integration_grid).
        ! As the grid is equidistant in t, the trapez sum
        !
        ! integral = (f_0/2 + f_N/2 + sum_i=1^{N-1} f_1)dt
        !
        ! has no f_N/2 term. Therefore, the right-endpoint (eta=eat_b, t=0) is
        ! excluded in eta_grid.
        integral = (sum(integrand*d_eta_dt) - 0.5_dp*integrand(1)*d_eta_dt(1))*dt
        deallocate (d_eta_dt)
    end function integrate_over_eta_grid

    function get_phi_integration_grid(fieldline) result(phi_grid)
        use utils, only: linspace
        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), allocatable :: phi_grid

        integer, parameter :: n_phi = 200
        real(dp), dimension(:), allocatable :: t
        real(dp) :: delta_phi, middle_phi

        allocate (phi_grid(n_phi), t(n_phi))
        ! transform integral to variable t -> phi = middle + delta * cos(t)
        call linspace(0.0_dp, pi, n_phi, t)
        delta_phi = 0.5_dp*(fieldline%phi_max(1) - fieldline%phi_max(2))
        middle_phi = 0.5_dp*(fieldline%phi_max(1) + fieldline%phi_max(2))
        phi_grid = middle_phi + delta_phi*cos(t)
        deallocate (t)
    end function get_phi_integration_grid

    subroutine allocate_eta_integrands(eta_grid, eta_integrand)
        real(dp), dimension(:), intent(in) :: eta_grid
        type(eta_integrand_t), intent(inout) :: eta_integrand

        integer :: n_eta

        if (allocated(eta_integrand%eta_grid)) deallocate (eta_integrand%eta_grid)
        if (allocated(eta_integrand%F1)) deallocate (eta_integrand%F1)
        if (allocated(eta_integrand%F2)) deallocate (eta_integrand%F2)
        if (allocated(eta_integrand%F3)) deallocate (eta_integrand%F3)
        if (allocated(eta_integrand%integral_lambda_over_B_squared)) then
            deallocate (eta_integrand%integral_lambda_over_B_squared)
        end if

        n_eta = size(eta_grid)
        allocate (eta_integrand%eta_grid, source=eta_grid)
        allocate (eta_integrand%F1(n_eta))
        allocate (eta_integrand%F2(n_eta))
        allocate (eta_integrand%F3(n_eta))
        allocate (eta_integrand%integral_lambda_over_B_squared(n_eta))

    end subroutine allocate_eta_integrands

    function cumint(x, defined_integral)
        real(dp), dimension(:), intent(in) :: x
        procedure(defined_integral_i) :: defined_integral
        real(dp), dimension(size(x)) :: cumint

        real(dp) :: x_start, x_end
        integer :: current, n_x

        n_x = size(x)
        cumint(1) = 0.0_dp ! as integral limits are same for first element
        do current = 2, n_x
            x_start = x(current - 1)
            x_end = x(current)
            cumint(current) = defined_integral(x_start, x_end) + cumint(current - 1)
        end do
    end function cumint

    function integral_one_over_lambda(phi_start, phi_end)
        use integrate, only: integrate_1d_substituted

        real(dp), intent(in) :: phi_start, phi_end
        real(dp) :: integral_one_over_lambda

        call integrate_1d_substituted(wrapper_one_over_lambda, &
                                      phi_start, &
                                      phi_end, &
                                      integral_one_over_lambda)

    end function integral_one_over_lambda

    function wrapper_one_over_lambda(phi)
        use fieldline_integrands, only: calc_lambda_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_one_over_lambda

        real(dp) :: theta, B, lambda

        theta = this_fieldline%get_theta(phi)
        call this_field%compute_B_mod(theta, phi, B)
        lambda = sqrt(calc_lambda_squared(B, this_eta))
        wrapper_one_over_lambda = 1.0_dp/lambda
    end function wrapper_one_over_lambda

    function wrapper_dBdtheta_over_B_squared_lambda(phi)
        use fieldline_integrands, only: calc_lambda_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_dBdtheta_over_B_squared_lambda

        real(dp) :: theta, B, dB_dx(3), dB_dtheta, lambda

        theta = this_fieldline%get_theta(phi)
        call this_field%compute_B_and_dB_dx(theta, phi, B, dB_dx)
        lambda = sqrt(calc_lambda_squared(B, this_eta))
        dB_dtheta = dB_dx(2)

        wrapper_dBdtheta_over_B_squared_lambda = dB_dtheta/(lambda*B**2.0_dp)
    end function wrapper_dBdtheta_over_B_squared_lambda

    function wrapper_lambda_dBdtheta_over_B_cubed(phi)
        use fieldline_integrands, only: calc_lambda_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_lambda_dBdtheta_over_B_cubed

        real(dp) :: theta, B, dB_dx(3), dB_dtheta, lambda

        theta = this_fieldline%get_theta(phi)
        call this_field%compute_B_and_dB_dx(theta, phi, B, dB_dx)
        lambda = sqrt(calc_lambda_squared(B, this_eta))
        dB_dtheta = dB_dx(2)

        wrapper_lambda_dBdtheta_over_B_cubed = lambda*dB_dtheta/(B**3.0_dp)
    end function wrapper_lambda_dBdtheta_over_B_cubed

    function wrapper_lambda_over_B_squared(phi)
        use fieldline_integrands, only: lambda_over_B_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_lambda_over_B_squared

        real(dp) :: theta

        theta = this_fieldline%get_theta(phi)
        wrapper_lambda_over_B_squared = lambda_over_B_squared(this_field, &
                                                              theta, &
                                                              phi, &
                                                              this_eta)
    end function wrapper_lambda_over_B_squared

end module shaing_callen_mod
