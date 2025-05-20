module fieldline_integrals
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    implicit none

contains

    subroutine calc_fieldline_integrals(field, fieldline)
        use integrate, only: integrate_1d
        type(fieldline_t), intent(inout) :: fieldline
        class(field_t), intent(in) :: field

        call integrate_1d(wrapper_lambda_over_B_squared, &
                          fieldline%phi_max(1), &
                          fieldline%phi_max(2), &
                          fieldline%I_hat)

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

    end subroutine calc_fieldline_integrals

    function get_theta(fieldline, phi) result(theta)
        type(fieldline_t), intent(in) :: fieldline
        real(dp) :: phi

        real(dp) :: theta

        theta = (phi - fieldline%phi_0)*fieldline%iota + fieldline%theta_0
    end function get_theta

end module fieldline_integrals
