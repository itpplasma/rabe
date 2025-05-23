module fieldline_integrals
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    implicit none

contains

    subroutine fourier_transform_over_label(field, fieldlines)
        use fourier, only: real_ft

        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        integer :: n_fieldlines

        integer :: current

        n_fieldlines = size(fieldlines)

        do current = 1, n_fieldlines
            call calc_fieldline_integrals(field, fieldlines(current))
        end do

    end subroutine fourier_transform_over_label

    subroutine calc_fieldline_integrals(field, fieldline)
        use integrate, only: integrate_1d
        type(fieldline_t), intent(inout) :: fieldline
        class(field_t), intent(in) :: field

        call integrate_1d(wrapper_lambda_over_B_squared, &
                          fieldline%phi_max(1), &
                          fieldline%phi_max(2), &
                          fieldline%well_average_lambda_b)

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
