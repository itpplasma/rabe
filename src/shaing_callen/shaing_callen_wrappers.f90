module shaing_callen_wrappers
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    implicit none

    class(field_t), allocatable :: this_field
    type(fieldline_t) :: this_fieldline
    type(fieldline_t) :: null_fieldline
    real(dp) :: this_eta
    real(dp), parameter :: null_eta = -1.0_dp

contains

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

    function wrapper_dBdtheta_over_lambda_cubed(phi)
        use fieldline_integrands, only: calc_lambda_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_dBdtheta_over_lambda_cubed

        real(dp) :: theta, B, dB_dx(3), dB_dtheta, lambda

        theta = this_fieldline%get_theta(phi)
        call this_field%compute_B_and_dB_dx(theta, phi, B, dB_dx)
        lambda = sqrt(calc_lambda_squared(B, this_eta))
        dB_dtheta = dB_dx(2)

        wrapper_dBdtheta_over_lambda_cubed = dB_dtheta/(lambda**3.0_dp)
    end function wrapper_dBdtheta_over_lambda_cubed

    function wrapper_dBdtheta_over_B_cubed(phi)
        real(dp), intent(in) :: phi
        real(dp) :: wrapper_dBdtheta_over_B_cubed

        real(dp) :: theta, B, dB_dx(3), dB_dtheta

        theta = this_fieldline%get_theta(phi)
        call this_field%compute_B_and_dB_dx(theta, phi, B, dB_dx)
        dB_dtheta = dB_dx(2)

        wrapper_dBdtheta_over_B_cubed = dB_dtheta/(B**3.0_dp)
    end function wrapper_dBdtheta_over_B_cubed

end module shaing_callen_wrappers
