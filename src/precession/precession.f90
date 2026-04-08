module precession
    use constants, only: dp, pi
    use field_base, only: field_t
    use field_base, only: field_3D_t
    use fieldline_mod, only: fieldline_t

    implicit none

    type :: integration_grid_t
        integer :: n_grid
        real(dp), dimension(:), allocatable :: eta
        real(dp), dimension(:), allocatable :: t
        real(dp), dimension(:), allocatable :: adiabatic_coef
        real(dp), dimension(:), allocatable :: bounce_coef
    end type integration_grid_t

    type, extends(fieldline_t) :: fieldline_with_minimum_t
        real(dp) :: phi_min
        real(dp) :: B_min
    end type fieldline_with_minimum_t

    type :: precession_t
        type(integration_grid_t) :: lower_grid
        type(integration_grid_t) :: upper_grid
        type(fieldline_with_minimum_t) :: fieldline
    end type precession_t

contains

 subroutine compute_precession_correction(field, fieldlines, l_c, Omega_hat, s_tor, correction)
        use field_instance, only: initialize_field_instance
        class(field_3D_t), intent(in) :: field
        class(fieldline_t), dimension(:) :: fieldlines
        real(dp), intent(in) :: l_c
        real(dp), intent(in) :: Omega_hat
        real(dp), intent(in) :: s_tor
        real(dp), intent(out) :: correction

        type(precession_t) :: precession
        real(dp) :: phi_bottom, B_bottom, lowest_B_max, eta_t, eta_c

        precession%fieldline%fieldline_t = get_fieldline_at_global_maximum(fieldlines)
        call find_magnetic_well_bottom(field, precession%fieldline, &
                                       phi_bottom, B_bottom)
        lowest_B_max = minval(precession%fieldline%B_max)
        eta_t = 1.0_dp/B_bottom
        eta_c = 1.0_dp/lowest_B_max
        precession%fieldline%phi_min = phi_bottom
        precession%fieldline%B_min = B_bottom
        call set_integration_grids(eta_t, eta_c, precession%lower_grid, &
                                   precession%upper_grid)
        call initialize_field_instance(field)
        call compute_bounce_integrals(field, precession%fieldline, s_tor, &
                                      precession%lower_grid)
        print *, "----------------------------------------------"
        print *, "Lower grid bounce integrals computed."
        print *, "----------------------------------------------"
        call compute_bounce_integrals(field, precession%fieldline, s_tor, &
                                      precession%upper_grid)

        correction = B_bottom

    end subroutine compute_precession_correction

    function get_fieldline_at_global_maximum(fieldlines) result(max_fieldline)
        class(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(fieldline_t) :: max_fieldline

        integer :: max_index

        max_index = maxloc(fieldlines%B_max(1), dim=1)
        max_fieldline = fieldlines(max_index)
    end function get_fieldline_at_global_maximum

    subroutine find_magnetic_well_bottom(field, fieldline, phi_bottom, B_bottom)
        use find_extrema, only: find_local_minima
        use, intrinsic :: ieee_arithmetic
        class(field_t), intent(in) :: field
        class(fieldline_t), intent(in) :: fieldline
        real(dp), intent(out) :: phi_bottom, B_bottom

        integer, parameter :: n_max = 10
        integer :: found_minima, of_smallest_B
        real(dp), dimension(n_max) :: phi_min, B_min

        call find_local_minima(B_along_fieldline, fieldline%phi_max, phi_min)

        found_minima = count(.not. ieee_is_nan(phi_min))
        if (found_minima < 1) then
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "error in find_magnetic_well_bottom: "
            print *, "Found no local minimum in provided interval!"
            print *, "theta_0: ", fieldline%theta_0
            print *, "phi_0: ", fieldline%phi_0
            print *, "interval: ", fieldline%phi_max
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            error stop
        elseif (found_minima > 2) then
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "warning in find_magnetic_well_bottom: "
            print *, "More than one local minimum in provided interval!"
            print *, "theta_0: ", fieldline%theta_0
            print *, "phi_0: ", fieldline%phi_0
            print *, "interval: ", fieldline%phi_max
            print *, "phi_min: ", phi_min(1:found_minima)
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
        end if

        call B_along_fieldline(phi_min(1:found_minima), B_min(1:found_minima))
        of_smallest_B = minloc(B_min(1:found_minima), dim=1)
        phi_bottom = phi_min(of_smallest_B)
        B_bottom = B_min(of_smallest_B)

    contains

        subroutine B_along_fieldline(phi, B)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: B

            real(dp), dimension(size(phi)) :: theta
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_mod(theta(idx), phi(idx), B(idx))
            end do
        end subroutine B_along_fieldline

    end subroutine find_magnetic_well_bottom

    subroutine set_integration_grids(eta_t, eta_c, lower_grid, upper_grid)
        use utils, only: linspace
        real(dp), intent(in) :: eta_t, eta_c
        type(integration_grid_t), intent(out) :: lower_grid, upper_grid

        real(dp) :: eta_mid
        real(dp) :: t_start, t_end
        integer, parameter :: n = 100
        real(dp), dimension(n) :: t

        if (allocated(lower_grid%eta)) deallocate (lower_grid%eta)
        if (allocated(lower_grid%t)) deallocate (lower_grid%t)

        allocate (lower_grid%eta(n))
        allocate (lower_grid%t(n))
        eta_mid = 0.5_dp*(eta_c + eta_t)
        t_start = 0.0_dp
        t_end = (eta_t - eta_mid)**0.25_dp
        call linspace(t_start, t_end, n, t)
        lower_grid%eta = eta_t - t**4.0_dp
        lower_grid%t = t
        lower_grid%n_grid = n

        if (allocated(upper_grid%eta)) deallocate (upper_grid%eta)
        if (allocated(upper_grid%t)) deallocate (upper_grid%t)

        allocate (upper_grid%eta(n))
        allocate (upper_grid%t(n))
        t_start = (eta_mid - eta_c)**0.25_dp
        t_end = 0.0_dp
        call linspace(t_start, t_end, n, t)
        upper_grid%eta = eta_c + t**4.0_dp
        upper_grid%t = t
        upper_grid%n_grid = n

    end subroutine set_integration_grids

    subroutine compute_bounce_integrals(field, fieldline, s_tor, grid)
        use field_instance, only: initialize_field_instance
        use params, only: params_init
        use bounce, only: trace_orbit_till_bounce
        use fieldline_integrands, only: calc_lambda_squared
        use constants, only: machine_eps
        class(field_3D_t), intent(in) :: field
        type(fieldline_with_minimum_t), intent(in) :: fieldline
        real(dp), intent(in) :: s_tor
        type(integration_grid_t), intent(inout) :: grid

        integer, parameter :: n_coef = 2
        integer :: idx, n_grid

        real(dp) :: theta_min, phi_min

        integer, parameter :: n_dim = 7
        real(dp), parameter :: cm2m = 1e-2_dp, gauss2tesla = 1e-4_dp
        real(dp) :: well_depth, I_j_estimate
        real(dp) :: deep_trapped_bounce_time, time_step
        real(dp), dimension(n_dim) :: z_template, z_start, z_end
        real(dp) :: lambda_squared
        real(dp) :: eta, t_weight
        real(dp) :: bounce_time_times_v_thermal
        real(dp) :: I_j
        real(dp) :: h_ctrvr(3), x(3)
        real(dp) :: dummy, dummy_vec(3)

        if (allocated(grid%adiabatic_coef)) deallocate (grid%adiabatic_coef)
        if (allocated(grid%bounce_coef)) deallocate (grid%bounce_coef)

        allocate (grid%adiabatic_coef(grid%n_grid))
        allocate (grid%bounce_coef(grid%n_grid))

        theta_min = fieldline%get_theta(fieldline%phi_min)
        phi_min = fieldline%phi_min
        x = [s_tor, theta_min, phi_min]
        call field%evaluate(x, dummy, dummy, dummy_vec, dummy_vec, h_ctrvr, dummy_vec)

        well_depth = 1.0_dp - fieldline%B_min*fieldline%eta_b
        I_j_estimate = 4.0_dp*sqrt(well_depth)/h_ctrvr(3)
        deep_trapped_bounce_time = 4.0_dp*pi/h_ctrvr(2)/sqrt(well_depth)
        time_step = deep_trapped_bounce_time/128.0_dp
        call params_init(fieldline%nfp, &
                         time_step, &
                         rlarm_in=0.0_dp)
        z_template(1) = s_tor
        z_template(2) = theta_min
        z_template(3) = phi_min
        z_template(4) = 1.0_dp
        z_template(5) = 0.0_dp
        z_template(6) = 0.0_dp
        z_template(7) = 0.0_dp

        do idx = 1, grid%n_grid
            t_weight = 4.0_dp*grid%t(idx)**3.0_dp
            if (abs(t_weight) < machine_eps) then
                grid%adiabatic_coef(idx) = 0.0_dp
                grid%bounce_coef(idx) = 0.0_dp
                cycle
            end if

            !> initialize orbit with certain eta
            z_start = z_template
            eta = grid%eta(idx)
            lambda_squared = calc_lambda_squared(fieldline%B_min, eta)
            if (abs(lambda_squared) < machine_eps) then
                print *, "idx: ", idx
                print *, "lambda_squared is zero!"
                print *, "B_min: ", fieldline%B_min
                print *, "eta: ", eta
                print *, "lambda_squared: ", lambda_squared
                error stop
            end if
            z_start(5) = sqrt(lambda_squared)

            !> initial trace to a turning point
            call trace_orbit_till_bounce(z_start, z_end)

            !> then reset time and integral and trace to until back at bounce
            !> i.e. one full bounce period
            z_start = z_end
            z_start(6:7) = 0.0_dp
            call trace_orbit_till_bounce(z_start, z_end)
            I_j = 0.5_dp*z_end(6)*cm2m/gauss2tesla ! I_j is only over half a bounce period
            bounce_time_times_v_thermal = z_end(7)*cm2m
            grid%adiabatic_coef(idx) = t_weight/(eta*I_j)
            grid%bounce_coef(idx) = t_weight*bounce_time_times_v_thermal
            print *, "idx: ", idx, "out of ", grid%n_grid
            print *, "Bounce time for eta = ", eta, ": ", bounce_time_times_v_thermal
            print *, "Deep trapped bounce time: ", deep_trapped_bounce_time
            print *, "Bounce integral I_j for eta = ", eta, ": ", I_j
            print *, "boundary layer estimate for I_j: ", I_j_estimate
        end do

    end subroutine compute_bounce_integrals

    function find_turning_points(field, fieldline, eta, reltol_in) result(phi_turning)
        use find_extrema, only: find_global_minimum
        class(field_t), intent(in) :: field
        type(fieldline_with_minimum_t), intent(in) :: fieldline
        real(dp), intent(in) :: eta
        real(dp), intent(in), optional :: reltol_in
        real(dp), dimension(2) :: phi_turning, extremum_locations
        real(dp) :: interval(2)

        real(dp) :: reltol

        if (present(reltol_in)) then
            reltol = reltol_in
        else
            reltol = 1e-4_dp
        end if

        ! turning points are shifted slightly away from the extrema to avoid
        ! numerical issues with the integrand
        interval(1) = fieldline%phi_max(1)
        interval(2) = fieldline%phi_min
        extremum_locations = find_global_minimum(lambda_squared_along_fieldline, interval, reltol)
        phi_turning(1) = extremum_locations(1) + (interval(2) - interval(1))*reltol

        interval(1) = fieldline%phi_min
        interval(2) = fieldline%phi_max(2)
        extremum_locations = find_global_minimum(lambda_squared_along_fieldline, interval, reltol)
        phi_turning(2) = extremum_locations(1) - (interval(2) - interval(1))*reltol

    contains

        subroutine lambda_squared_along_fieldline(phi, lambda_squared)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: lambda_squared

            real(dp), dimension(size(phi)) :: theta, B
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_mod(theta(idx), phi(idx), B(idx))
            end do
            lambda_squared = abs(1.0_dp - B*eta)
        end subroutine lambda_squared_along_fieldline

    end function find_turning_points

end module precession
