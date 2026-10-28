module precession_reference
    use constants, only: dp
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t
    use grid_mod, only: integration_grid_t
    use grid_mod, only: fieldline_for_precession_t
    use interpolate, only: evaluate_splines_1d, evaluate_splines_1d_many

    implicit none

    private
    class(integration_grid_t), allocatable :: grid
    logical :: grid_initialized = .false.

    public :: integrate_radial_drift
    public :: get_fieldline_at_global_maximum
    public :: find_turning_points
    public :: evaluate_grid_splines

contains

    subroutine integrate_radial_drift(grid, fieldline, integral)
        use odeint_allroutines_sub, only: odeint_allroutines
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
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), intent(out) :: dydt

        real(dp) :: radial_drift
        call get_radial_drift(t, dydt(1))

    end subroutine rhs

    subroutine initialize_grid_instance(grid_in)
        class(integration_grid_t), intent(in) :: grid_in
        if (grid_initialized) then
            deallocate (grid)
        end if
        allocate (grid, source=grid_in)
        grid_initialized = .true.
    end subroutine initialize_grid_instance

    subroutine get_radial_drift(t, radial_drift)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: radial_drift

        if (.not. grid_initialized) then
        error stop "Error: get_radial_drift called before grid instance is initialized."
        end if
        call evaluate_splines_1d(grid%radial_drift_weighted_spline, t, radial_drift)
    end subroutine get_radial_drift

    subroutine get_bounce_time(t, bounce_time)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: bounce_time

        if (.not. grid_initialized) then
         error stop "Error: get_bounce_time called before grid instance is initialized."
        end if
        call evaluate_splines_1d(grid%bounce_time_weighted_spline, t, bounce_time)
    end subroutine get_bounce_time

    subroutine get_I_j(t, I_j)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: I_j

        if (.not. grid_initialized) then
            error stop "Error: get_I_j called before grid instance is initialized."
        end if
        call evaluate_splines_1d(grid%I_j_spline, t, I_j)
    end subroutine get_I_j

    subroutine get_poloidal_drift(t, poloidal_drift)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: poloidal_drift

        if (.not. grid_initialized) then
            print *, "Error: get_poloidal_drift called before ", &
                "grid instance is initialized."
            error stop
        end if
        call evaluate_splines_1d(grid%poloidal_drift_weighted_spline, t, poloidal_drift)
    end subroutine get_poloidal_drift

    function get_fieldline_at_global_maximum(fieldlines) result(max_fieldline)
        class(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(fieldline_t) :: max_fieldline

        integer :: max_index

        max_index = maxloc(fieldlines%B_max(1), dim=1)
        max_fieldline = fieldlines(max_index)
    end function get_fieldline_at_global_maximum

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

    subroutine evaluate_grid_splines(grid, t_eval, I_j_eval, bounce_coef_eval, &
                                     radial_drift_eval, poloidal_drift_eval)
        type(integration_grid_t), intent(in) :: grid
        real(dp), intent(in) :: t_eval(:)
        real(dp), intent(out) :: I_j_eval(:)
        real(dp), intent(out) :: bounce_coef_eval(:)
        real(dp), intent(out), optional :: radial_drift_eval(:)
        real(dp), intent(out), optional :: poloidal_drift_eval(:)

        call evaluate_splines_1d_many(grid%I_j_spline, t_eval, I_j_eval)
        call evaluate_splines_1d_many(grid%bounce_time_weighted_spline, t_eval, bounce_coef_eval)
        if (present(radial_drift_eval)) then
            call evaluate_splines_1d_many(grid%radial_drift_weighted_spline, t_eval, radial_drift_eval)
        end if
        if (present(poloidal_drift_eval)) then
            call evaluate_splines_1d_many(grid%poloidal_drift_weighted_spline, t_eval, poloidal_drift_eval)
        end if
    end subroutine evaluate_grid_splines

end module precession_reference
