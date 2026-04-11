module precession
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use constants, only: dp, pi
    use field_base, only: field_t
    use field_base, only: field_3D_t
    use fieldline_mod, only: fieldline_t
    use grid_mod, only: integration_grid_t
    use grid_mod, only: fieldline_for_precession_t
    use grid_mod, only: set_integration_grids
    use grid_mod, only: compute_bounce_integrals
    use grid_mod, only: set_splines
    use interpolate, only: SplineData1D, construct_splines_1d, evaluate_splines_1d_many
    use fieldline_integrals, only: modes_t

    implicit none

contains

    subroutine compute_precession_correction(field, fieldlines_in, l_c, Omega_hat, s_tor, correction)
        use field_instance, only: initialize_field_instance
        use fourier, only: real_ft
        use surface_average_mod, only: calc_surface_averages
        use surface_average_mod, only: surface_average_t
        use splines_instance, only: initialize_splines
        use splines_instance, only: initialize_prefactor
        use splines_instance, only: initialize_radial_drift_spline
        use splines_instance, only: get_flux_mode
        class(field_3D_t), intent(in) :: field
        class(fieldline_t), dimension(:) :: fieldlines_in
        real(dp), intent(in) :: l_c
        real(dp), intent(in) :: Omega_hat
        real(dp), intent(in) :: s_tor
        real(dp), intent(out) :: correction

        type(fieldline_for_precession_t), dimension(:), allocatable :: fieldlines
        real(dp) :: phi_bottom, B_bottom, lowest_B_max, eta_t, eta_c
        type(integration_grid_t) :: grid

        real(dp), dimension(:, :), allocatable :: radial_drift_weighted
        real(dp), dimension(:, :), allocatable :: radial_drift_cos
        real(dp), dimension(:, :), allocatable :: radial_drift_sin
        real(dp), dimension(:), allocatable :: bounce_time_weighted
        real(dp), dimension(:), allocatable :: I_j
        real(dp), dimension(:), allocatable :: poloidal_drift_weighted
        real(dp), dimension(:), allocatable :: electric_drift_weighted

        real(dp), dimension(:), allocatable :: flux_mode
        real(dp) :: M_pol, N_tor, nfp, iota
        real(dp) :: mode_factor
        type(surface_average_t) :: average

        integer :: n_fieldlines
        integer :: n_modes
        integer :: idx

        n_fieldlines = size(fieldlines_in)
        allocate (fieldlines(n_fieldlines))
        fieldlines%fieldline_t = fieldlines_in

        eta_c = 1.0_dp/get_smallest_maximum(fieldlines)
        call set_fieldline_minima(field, fieldlines)
        eta_t = 1.0_dp/get_biggest_minimum(fieldlines)
        call set_integration_grids(eta_t, eta_c, grid)
        fieldlines%grid = grid

        call initialize_field_instance(field)
        allocate (radial_drift_weighted(n_fieldlines, grid%n))
        allocate (bounce_time_weighted(grid%n))
        allocate (I_j(grid%n))
        allocate (poloidal_drift_weighted(grid%n))
        do idx = 1, n_fieldlines
            call compute_bounce_integrals(field, &
                                          fieldlines(idx), &
                                          s_tor, &
                                          fieldlines(idx)%grid)
            radial_drift_weighted(idx, :) = fieldlines(idx)%grid%radial_drift_weighted
 bounce_time_weighted = bounce_time_weighted + fieldlines(idx)%grid%bounce_time_weighted
            I_j = I_j + fieldlines(idx)%grid%I_j
            poloidal_drift_weighted = poloidal_drift_weighted + fieldlines(idx)%grid%poloidal_drift_weighted
        end do
        electric_drift_weighted = Omega_hat*bounce_time_weighted

        n_modes = n_fieldlines/2 + 1
        allocate (radial_drift_cos(n_fieldlines, grid%n))
        allocate (radial_drift_sin(n_fieldlines, grid%n))
        do idx = 1, grid%n
            call real_ft(fieldlines%xi_0, &
                         radial_drift_weighted(:, idx), &
                         radial_drift_cos(:, idx), &
                         radial_drift_sin(:, idx))
        end do

        call initialize_splines(grid%t, &
                                grid%eta, &
                                I_j, &
                                electric_drift_weighted, &
                                poloidal_drift_weighted)

        allocate (flux_mode(n_modes))
        flux_mode(1) = 0.0_dp
        do idx = 2, n_modes
            mode_factor = real(idx - 1, dp)*l_c*nfp/(M_pol*iota - N_tor)
            call initialize_prefactor(mode_factor)
            call initialize_radial_drift_spline(grid%t, radial_drift_sin(idx, :))
            call get_flux_mode(flux_mode(idx))
        end do

        call calc_surface_averages(fieldlines, average)
        correction = pi*sum(flux_mode)/average%normalization

    end subroutine compute_precession_correction

    function get_smallest_maximum(fieldlines) result(smallest_maximum)
        class(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: smallest_maximum

        real(dp), dimension(size(fieldlines)) :: smallest_maximum_per_line
        integer :: idx

        do idx = 1, size(fieldlines)
            smallest_maximum_per_line(idx) = minval(fieldlines(idx)%B_max)
        end do
        smallest_maximum = minval(smallest_maximum_per_line)

        if (smallest_maximum <= 0.0_dp) then
            error stop "Error: Smallest maximum is not positiv."
        end if
        if (ieee_is_nan(smallest_maximum)) then
            error stop "Error: Smallest maximum is NaN."
        end if

    end function get_smallest_maximum

    subroutine set_fieldline_minima(field, fieldlines)
        class(field_t), intent(in) :: field
        class(fieldline_for_precession_t), dimension(:), intent(inout) :: fieldlines

        integer :: idx

        do idx = 1, size(fieldlines)
            call find_magnetic_well_bottom(field, &
                                           fieldlines(idx)%fieldline_t, &
                                           fieldlines(idx)%phi_min, &
                                           fieldlines(idx)%B_min)
        end do
    end subroutine set_fieldline_minima

    function get_biggest_minimum(fieldlines) result(biggest_minimum)
        class(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: biggest_minimum

        biggest_minimum = maxval(fieldlines%B_min)

        if (biggest_minimum <= 0.0_dp) then
            error stop "Error: Biggest minimum is not positiv."
        end if
        if (ieee_is_nan(biggest_minimum)) then
            error stop "Error: Biggest minimum is NaN."
        end if
    end function get_biggest_minimum

    subroutine integrate_radial_drift(grid, fieldline, integral)
        use odeint_allroutines_sub, only: odeint_allroutines
        use grid_instance, only: initialize_grid_instance
        type(integration_grid_t), intent(in) :: grid
        type(fieldline_t), intent(in) :: fieldline
        real(dp), intent(out) :: integral

        integer, parameter :: ndim = 1
        real(dp) :: t_start, t_end
        real(dp), parameter :: relerr = 1e-8_dp
        real(dp), dimension(ndim) :: y

        t_start = grid%t(1)
        t_end = grid%t(size(grid%t))

        call initialize_grid_instance(grid)
        y(1) = 0.0_dp
        call odeint_allroutines(y, ndim, t_start, t_end, relerr, rhs)

        integral = y(1)
    end subroutine integrate_radial_drift

    subroutine rhs(t, y, dydt)
        use grid_instance, only: get_radial_drift
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), intent(out) :: dydt

        real(dp) :: radial_drift
        call get_radial_drift(t, dydt(1))

    end subroutine rhs

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

    function find_turning_points(field, fieldline, eta, reltol_in) result(phi_turning)
        use find_extrema, only: find_global_minimum
        class(field_t), intent(in) :: field
        type(fieldline_for_precession_t), intent(in) :: fieldline
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
