
module splines_instance
    use constants, only: dp
    use interpolate, only: evaluate_splines_1d
    use interpolate, only: SplineData1D
    use interpolate, only: construct_splines_1d

    implicit none

    private
    type(SplineData1D) :: I_j_spl
    type(SplineData1D) :: magnetic_drift_spl
    type(SplineData1D) :: electric_drift_spl
    type(SplineData1D) :: radial_drift_mode_spl
    type(SplineData1D) :: eta_spl
    real(dp) :: prefactor

    logical :: splines_initialized = .false.
    logical, parameter :: periodic = .false.
    integer, parameter :: order = 3

    public :: get_flux_mode
    public :: get_radial_drift_mode
    public :: get_magnetic_drift
    public :: get_I_j
    public :: get_electric_drift
    public :: initialize_splines
    public :: initialize_radial_drift_spline
    public :: initialize_prefactor

contains
    subroutine get_flux_mode(flux_mode)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        real(dp), intent(out) :: flux_mode
        if (.not. splines_initialized) then
          print *, "Error: get_flux_mode called before splines instance is initialized."
            error stop
        end if
        flux_mode = ieee_value(1.0_dp, ieee_quiet_nan)
    end subroutine get_flux_mode

    subroutine get_radial_drift_mode(t, radial_drift)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: radial_drift

        if (.not. splines_initialized) then
     error stop "Error: get_radial_drift called before splines instance is initialized."
        end if
        call evaluate_splines_1d(radial_drift_mode_spl, t, radial_drift)
    end subroutine get_radial_drift_mode

    subroutine get_magnetic_drift(t, magnetic_drift)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: magnetic_drift

        if (.not. splines_initialized) then
   error stop "Error: get_magnetic_drift called before splines instance is initialized."
        end if
        call evaluate_splines_1d(magnetic_drift_spl, t, magnetic_drift)
    end subroutine get_magnetic_drift

    subroutine get_I_j(t, I_j)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: I_j

        if (.not. splines_initialized) then
            error stop "Error: get_I_j called before splines instance is initialized."
        end if
        call evaluate_splines_1d(I_j_spl, t, I_j)
    end subroutine get_I_j

    subroutine get_electric_drift(t, electric_drift)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: electric_drift

        if (.not. splines_initialized) then
            print *, "Error: get_electric_drift called before ", &
                "splines instance is initialized."
            error stop
        end if
        call evaluate_splines_1d(electric_drift_spl, t, electric_drift)
    end subroutine get_electric_drift

    subroutine initialize_splines(t, eta, I_j, magnetic_drift, electric_drift)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        real(dp), dimension(:), intent(in) :: t
        real(dp), dimension(:), intent(in) :: eta
        real(dp), dimension(:), intent(in) :: I_j
        real(dp), dimension(:), intent(in) :: magnetic_drift
        real(dp), dimension(:), intent(in) :: electric_drift

        integer :: start, n

        if (splines_initialized) then
            error stop "Error: initialize_splines_instance called more than once."
        end if

        n = size(t)
        if (ieee_is_nan(I_j(1))) then
            start = 2
        else
            start = 1
        end if
        call construct_splines_1d(t(start), t(n), I_j(start:n), &
                                  order, periodic, I_j_spl)

        start = 1
        call construct_splines_1d(t(start), t(n), magnetic_drift(start:n), &
                                  order, periodic, magnetic_drift_spl)
        call construct_splines_1d(t(start), t(n), electric_drift(start:n), &
                                  order, periodic, electric_drift_spl)
        call construct_splines_1d(t(start), t(n), eta(start:n), &
                                  order, periodic, eta_spl)

        splines_initialized = .true.
    end subroutine initialize_splines

    subroutine initialize_radial_drift_spline(t, radial_drift_mode)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        real(dp), dimension(:), intent(in) :: t
        real(dp), dimension(:), intent(in) :: radial_drift_mode

        integer :: start, n

        n = size(t)
        start = 1
        call construct_splines_1d(t(start), t(n), radial_drift_mode(start:n), &
                                  order, periodic, radial_drift_mode_spl)
    end subroutine initialize_radial_drift_spline

    subroutine initialize_prefactor(prefactor_in)
        real(dp), intent(in) :: prefactor_in
        prefactor = prefactor_in
    end subroutine initialize_prefactor
end module splines_instance
