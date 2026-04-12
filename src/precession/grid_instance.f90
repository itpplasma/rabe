
module grid_instance
    use constants, only: dp
    use field_base, only: field_3D_t
    use grid_mod, only: integration_grid_t
    use grid_mod, only: evaluate_grid_splines
    use interpolate, only: evaluate_splines_1d

    implicit none

    private
    class(integration_grid_t), allocatable :: grid
    logical :: grid_initialized = .false.

    public :: initialize_grid_instance, get_radial_drift

contains

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

end module grid_instance
