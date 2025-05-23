module fieldline_integrals
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    implicit none

    type :: modes_t
        real(dp), dimension(:), allocatable :: cos_coeffs, sin_coeffs
        real(dp), dimension(:), allocatable :: mode_number
    end type modes_t

contains

    subroutine fourier_transform_over_label(field, fieldlines, radial_drift_modes)
        use fourier, only: real_ft

        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        type(modes_t), intent(out) :: radial_drift_modes

        integer :: n_fieldlines, n_modes
        integer :: current
        real(dp), dimension(size(fieldlines)) :: label, radial_drift

        n_fieldlines = size(fieldlines)
        n_modes = n_fieldlines/2 + 1

        do current = 1, n_fieldlines
            call calc_fieldline_integrals(field, fieldlines(current))
        end do

        label = fieldlines(:)%theta_0
        radial_drift = fieldlines(:)%well_average_radial_drift_velocity

        allocate (radial_drift_modes%cos_coeffs(n_modes))
        allocate (radial_drift_modes%sin_coeffs(n_modes))
        allocate (radial_drift_modes%mode_number(n_modes))

        call real_ft(label, &
                     radial_drift, &
                     radial_drift_modes%cos_coeffs, &
                     radial_drift_modes%sin_coeffs)

    end subroutine fourier_transform_over_label

    subroutine calc_fieldline_integrals(field, fieldline)
        use integrate, only: integrate_1d
        type(fieldline_t), intent(inout) :: fieldline
        class(field_t), intent(in) :: field

        call integrate_1d(wrapper_radial_drift_velocity, &
                          fieldline%phi_max(1), &
                          fieldline%phi_max(2), &
                          fieldline%well_average_radial_drift_velocity)

    contains

        function wrapper_lambda_over_B_squared(phi)
            use fieldline_integrands, only: lambda_over_B_squared

            real(dp), intent(in) :: phi
            real(dp) :: wrapper_lambda_over_B_squared

            real(dp) :: theta

            theta = get_theta(fieldline, phi)
            wrapper_lambda_over_B_squared = lambda_over_B_squared(field, &
                                                                  theta, &
                                                                  phi, &
                                                                  fieldline%eta_b)
        end function wrapper_lambda_over_B_squared

        function wrapper_radial_drift_velocity(phi)
            use fieldline_integrands, only: radial_drift_velocity

            real(dp), intent(in) :: phi
            real(dp) :: wrapper_radial_drift_velocity

            real(dp) :: theta

            theta = get_theta(fieldline, phi)
            wrapper_radial_drift_velocity = radial_drift_velocity(field, &
                                                                  theta, &
                                                                  phi, &
                                                                  fieldline%eta_b)
        end function wrapper_radial_drift_velocity

    end subroutine calc_fieldline_integrals

    function get_theta(fieldline, phi) result(theta)
        type(fieldline_t), intent(in) :: fieldline
        real(dp) :: phi

        real(dp) :: theta

        theta = (phi - fieldline%phi_0)*fieldline%iota + fieldline%theta_0
    end function get_theta

end module fieldline_integrals
